library(RSNNS)

data(snnsData)

inputs <- snnsData$eight_016.pat[,inputColumns(snnsData$eight_016.pat)]
outputs <- snnsData$eight_016.pat[,outputColumns(snnsData$eight_016.pat)]

par(mfrow=c(1,2))

modelElman <- elman(inputs, outputs, size=8, learnFuncParams=c(0.1), maxit=1000)
modelElman

modelJordan <- jordan(inputs, outputs, size=8, learnFuncParams=c(0.1), maxit=1000)
modelJordan

plotIterativeError(modelElman)
plotIterativeError(modelJordan)

summary(modelElman)
summary(modelJordan)