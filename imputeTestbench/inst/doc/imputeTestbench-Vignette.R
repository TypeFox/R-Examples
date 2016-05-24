## ------------------------------------------------------------------------
datax <- c(1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5)

## ------------------------------------------------------------------------
library(imputeTestbench)

## ------------------------------------------------------------------------
q <- impute_errors(datax)
q
plot_errors(q)

## ------------------------------------------------------------------------
#aa <- append_method(existing_method = q,dataIn= datax,missPercentFrom = 10, missPercentTo = 80, interval = 10, MethodPath = "source('~/imputeTestbench/R/inter.R')", MethodName = "Random")

#aa
#plot_errors(aa)

## ------------------------------------------------------------------------
#bb <- append_method(existing_method = aa, dataIn= datax,missPercentFrom = 10, missPercentTo = 80, interval = 10, MethodPath = "source('~/imputeTestbench/R/PSFimpute.R')", MethodName = "PSFimpute")

#bb
#plot_errors(bb)

