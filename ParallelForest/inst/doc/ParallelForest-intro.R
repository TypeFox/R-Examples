## ------------------------------------------------------------------------
library(ParallelForest)

data(low_high_earners)       # cleaned and prepared training dataset
data(low_high_earners_test)  # cleaned and prepared testing dataset

## ------------------------------------------------------------------------
fforest = grow.forest(Y~., data=low_high_earners)

## ------------------------------------------------------------------------
fforest["min_node_obs"]

## ------------------------------------------------------------------------
fforest["max_depth"]

## ------------------------------------------------------------------------
fforest["numsamps"]

## ------------------------------------------------------------------------
fforest["numvars"]

## ------------------------------------------------------------------------
fforest["numboots"]

## ------------------------------------------------------------------------
fforest_samepred = predict(fforest, low_high_earners)
pctcorrect_samepred = sum(low_high_earners$Y==fforest_samepred)/nrow(low_high_earners)
print(pctcorrect_samepred)

## ------------------------------------------------------------------------
fforest_newpred = predict(fforest, low_high_earners_test)
pctcorrect_newpred = sum(low_high_earners$Y==fforest_newpred)/nrow(low_high_earners)
print(pctcorrect_newpred)

