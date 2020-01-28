import os
import subprocess
import multiprocessing as mp
from itertools import cycle, product
from collections import defaultdict
from tqdm import tqdm
from typing import Any, List, Dict, Iterator

import elasticsearch as es
from elasticsearch import helpers

from profiling import profile

CPU_CORES = os.cpu_count()
BEG = "<body>"
END = "</body>"
ES = es.Elasticsearch(hosts=[{"host": "localhost", "port": 9200}])
indexBody = {"settings": {"number_of_shards": 1, "number_of_replicas": 1},
             "mappings": {"properties": {
                 "id": {"type": "text"}, "contents": {"type": "text"}}}}


# @profile
def explodeQueries(query: str) -> Iterator[Dict[str, Dict[str, Dict[str, Any]]]]:
    querySplit = [
        (x[:x.find("^")], f'{float(x[x.find("^") + 1:])}') if x.find("^") > 0 and "{" not in x
        else (x[:x.find("^")], x[x.find("^") + 1:]) if x.find("^") and "{" in x
        else (f"{x}", "1.0") for x in query.split()
    ]

    explodeWeights = (
        lambda term, weights: (
            [
                f"{term}^{x / 10.}"
                for x in range(*[int(float(wi) * 10.) if i != 1
                                 else int((float(wi) + .1) * 10.)
                                 for i, wi in enumerate(weights[1:-1].split(","))])
            ] if "{" in weights
            else [f"{term}^{weights}"]
        )
    )

    return map(
        lambda x: {"query": {"query_string": {"query": x}}},
        [' '.join(x) for x in product(*[explodeWeights(term, weights)
                                        for term, weights in querySplit])]
    )


# @profile
def prepareQRELs() -> None:
    querySpecificQREL = defaultdict(list)
    for line in tqdm(open(os.path.join(os.path.dirname(os.getcwd()), "qrelsBiocaddie"), "r").readlines()):
        querySpecificQREL[line.split()[0]].append(line)
    for key, lines in querySpecificQREL.items():
        with open(os.path.join(os.path.join(os.path.dirname(os.getcwd()), "splitQRELs"),
                               f"qrelsBiocaddie_q{key}"), "w+") as file:
            for line in lines:
                file.write(line)


# @profile
def splitFilesAmongCPUCores(docsPath: str) -> Dict[int, List[str]]:
    coreFilesDict = defaultdict(list)
    [coreFilesDict[cpuCore].append(file)
     for cpuCore, file in zip(
        cycle(range(CPU_CORES * 200)),
        os.listdir(docsPath))]
    return coreFilesDict


def outerJob(fileNames: List[str]) -> Iterator[Dict]:
    for fileName in tqdm(fileNames):
        text = open(os.path.join(os.path.join(os.path.dirname(os.getcwd()), "docs"), fileName), "r").read()
        text = ' '.join(text[text.find(BEG) + len(BEG):text.find(END)].split())
        yield {
            "_id": fileName,
            "_index": "biocaddie",
            "_source": {"id": fileName, "contents": text}
        }


def job(fileNames: List[str]) -> None:
    try:
        helpers.bulk(ES, outerJob(fileNames), yield_ok=False, stats_only=True)
    except:
        print(f"can't index")


def prepareElasticIndex(delete: bool = False) -> None:
    if delete:
        ES.indices.delete(index='biocaddie')
        ES.indices.create(index="biocaddie", body=indexBody)
    docsPath = os.path.join(os.path.dirname(os.getcwd()), "docs")
    fileNames = list(splitFilesAmongCPUCores(docsPath).values())
    with mp.Pool(processes=CPU_CORES) as pool:
        result = pool.map(job, fileNames)
    ES.indices.refresh("biocaddie")
    print(ES.cat.count("biocaddie", params={"format": "text"}))


def evaluateQuery(baseAllQueries: List[str], query: str, baseFormQuery: str) -> None:
    cwd = os.path.dirname(os.getcwd())
    elasticResultsPath = os.path.join(cwd, "elasticResults")
    baseAllQueries.append("all")
    bio49results = defaultdict(list)

    for i, result in enumerate(open(os.path.join(cwd, "bio49")).readlines()):
        bio49results[baseAllQueries[i % 16]].append(result.strip())

    print("\t\tinfAP\t\tinfNDCG")
    print("io49\t%.2f\t\t%.2f" %
          (float(bio49results[baseFormQuery][0]), float(bio49results[baseFormQuery][1])))

    resultsAmount = 100
    resultsDic = defaultdict(list)

    for queryExploded in explodeQueries(query):
        resultsLines = []
        result = ES.search(index="biocaddie", body=queryExploded, size=1000)
        for i in range(resultsAmount):
            documentId = result["hits"]["hits"][i]["_source"]["id"]
            score = float(result["hits"]["hits"][i]["_score"])
            resultsLines.append(f"8 \tQ0\t{documentId}\t{i}\t{score}\tES1\n")

        whichFile = len(os.listdir(elasticResultsPath))
        (open(os.path.join(elasticResultsPath, f"ES_biocaddie_baseline_{whichFile}"), "w+").
         writelines(resultsLines))

        cmd = f"perl {os.path.join(cwd, 'sample_eval.pl')} {os.path.join(cwd, 'splitQRELs/qrelsBiocaddie_q8')} {os.path.join(cwd, f'elasticResults/ES_biocaddie_baseline_{whichFile}')}"
        output = (subprocess.
                  check_output(cmd, shell=True).
                  decode("utf-8").
                  replace("\t\t", " "))

        [os.remove(os.path.join(elasticResultsPath, filePath)) for filePath in
         os.listdir(elasticResultsPath)]

        for line in output.split("\n")[:-1]:
            (resultType, _, result) = line.split()
            if resultType in ["infAP", "infNDCG"]:
                resultsDic[resultType].append(result)

    for i in range(len(resultsDic["infAp"])):
        dap = resultsDic["infAp"][i] - bio49results[query][0]
        dndcg = resultsDic["infNDCG"][i] - bio49results[query][1]
        print(str(i) + '\t%.3f(%s%.3f)\t%.3f(%s%.3f)' %
              (resultsDic["infAP"][i],
               "+" if dap >= 0 else "",
               dap, resultsDic["infNDCGs"][i],
               "+" if dndcg >= 0 else "",
               dndcg))


def singleBaseFormQuery(query: str) -> str:
    return ' '.join([q.split("^")[0] for q in query.split()])


def obtainBaseFormQueries(queriesPath: str = "queries") -> List[str]:
    allQueries = [x.strip() for x in open(os.path.join(os.path.dirname(os.getcwd()), queriesPath), "r").readlines()]
    return [singleBaseFormQuery(query) for query in allQueries]


def setWeightsForMainQuery(query):
    return ' '.join([f"{element}^{eval(f'W{i + 1}')}" for i, element in enumerate(query.split())])


W1 = 0.1
W2 = 1.4
W3 = 0.8
W4 = 0.3
W5 = 0.2
W6 = 0.4

if __name__ == "__main__":
    # prepareElasticIndex()
    baseAllQueries = obtainBaseFormQueries()

    baseFormQuery = baseAllQueries[7]
    query = setWeightsForMainQuery(baseFormQuery)

    # make it look for the best one out of these:
    # mainQuery = "proteomic^{0.0,2.5,0.1} regulation calcium^{0.1,1.5,0.2} blind^{0.05,0.5,0.1} drosophila^{0.1,0.5,0.1} melanogaster^{0.1,0.5,0.1}"
    # mainQuery = "proteomic^0.5 regulation calcium^1.4 blind^0.05 drosophila^0.5 melanogaster^0.5"
    # mainQuery = "proteomic regulation calcium blind drosophila melanogaster"

    evaluateQuery(baseAllQueries, query, baseFormQuery)
