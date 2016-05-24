library(RSNNS)

data(snnsData)

inputs <- snnsData$laser_1000.pat[,inputColumns(snnsData$laser_1000.pat)]
outputs <- snnsData$laser_1000.pat[,outputColumns(snnsData$laser_1000.pat)]

patterns <- splitForTrainingAndTest(inputs, outputs, ratio=0.15)

model <- elman(patterns$inputsTrain, patterns$targetsTrain, size=c(8,8), learnFuncParams=c(0.1), maxit=500,
    inputsTest=patterns$inputsTest, targetsTest=patterns$targetsTest, linOut=FALSE)

modelJordan <- jordan(patterns$inputsTrain, patterns$targetsTrain, size=c(8), learnFuncParams=c(0.1), maxit=100,
                inputsTest=patterns$inputsTest, targetsTest=patterns$targetsTest, linOut=FALSE)

#modelMlp <- mlp(patterns$inputsTrain, patterns$targetsTrain, initFuncParams=c(-0.3,0.3),size=c(8), learnFuncParams=c(0.05), maxit=500,
#                inputsTest=patterns$inputsTest, targetsTest=patterns$targetsTest, linOut=TRUE)

names(model)
#model$IterativeFitError
#model$fitted.values
#model$fittedTestValues

par(mfrow=c(3,3))

plotIterativeError(model)
plotIterativeError(modelJordan)
#plotIterativeError(modelMlp)

plotRegressionError(patterns$targetsTrain, model$fitted.values, main="Regression Plot Fit")
plotRegressionError(patterns$targetsTest, model$fittedTestValues, main="Regression Plot Test")
hist(model$fitted.values - patterns$targetsTrain, col="lightblue", main="Error Histogram Fit")

#model$IterativeFitError[length(model$IterativeFitError)]

plot(inputs, type="l")

plot(inputs[1:100], type="l")
lines(outputs[1:100], col="red")
lines(model$fitted.values[1:100], col="green")