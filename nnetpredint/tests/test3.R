# Example 3: Using the rsnns object trained by RSNNS package
library(RSNNS)
data(iris)
#shuffle the vector
iris <- iris[sample(1:nrow(iris),length(1:nrow(iris))),1:ncol(iris)]
irisValues <- iris[,1:4]
irisTargets <- decodeClassLabels(iris[,5])[,'setosa']

iris <- splitForTrainingAndTest(irisValues, irisTargets, ratio=0.15)
iris <- normTrainingAndTestSet(iris)
model <- mlp(iris$inputsTrain, iris$targetsTrain, size=5, learnFuncParams=c(0.1),
	maxit=50, inputsTest=iris$inputsTest, targetsTest=iris$targetsTest)
predictions <- predict(model,iris$inputsTest)


# Generating prediction intervals
library(nnetpredint)

# S3 Method for rsnns class prediction intervals
xTrain <- iris$inputsTrain
yTrain <- iris$targetsTrain
newData <- iris$inputsTest
yPredInt <- nnetPredInt(model, xTrain, yTrain, newData)
print(yPredInt[1:20,])

