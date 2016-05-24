## ------------------------------------------------------------------------
datax <- c(1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5)

## ------------------------------------------------------------------------
library(PredictTestbench)

## ------------------------------------------------------------------------
q <- prediction_errors(datax, nextVal = 10)
q
plot_predictions(q)

## ------------------------------------------------------------------------
#aa <- prediction_errors(dataIn= datax, nextVal = 10, MethodPath = "source('~/PredictTestbench/R/inter.R')", MethodName = "Proposed_Method")

#aa
#plot_predictions(aa)

## ------------------------------------------------------------------------
#bb <- prediction_append(existing_method = aa, dataIn= datax, nextVal = 10, MethodPath = "source('~/imputeTestbench/R/PSFimpute.R')", MethodName = "PSFimpute")

#bb
#plot_predictions(bb)

## ------------------------------------------------------------------------
# cc <- prediction_remove (existing_method = bb, index_number = 1)
# cc
# plot_predictions(cc)

