library(RSNNS)

set.seed(2)

inputs <- as.matrix(seq(0,10,0.1))
outputs <- as.matrix(sin(inputs) + runif(inputs*0.2))
#outputs <- as.matrix(sin(inputs))
outputs <- normalizeData(outputs, "0_1")

model <- rbf(inputs, outputs, size=40, maxit=1000, initFuncParams=c(0, 1, 0, 0.01, 0.01), 
    learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)

par(mfrow=c(2,1))

plotIterativeError(model)

plot(inputs, outputs)
lines(inputs, fitted(model), col="green")
