guild run query-evaluator \
W1=loguniform[0.0:2.0] \
W2=loguniform[0.0:2.0] \
W3=loguniform[0.0:2.0] \
W4=loguniform[0.0:2.0] \
W5=loguniform[0.0:2.0] \
W6=loguniform[0.0:2.0] \
--optimizer gp \
--max-trials 100 \
--label trec-test