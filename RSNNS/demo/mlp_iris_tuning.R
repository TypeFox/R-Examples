library(RSNNS)

set.seed(2)

data(iris)

#shuffle the vector
iris <- iris[sample(nrow(iris)),]

irisValues <- iris[,1:4]
#irisTargets <- decodeClassLabels(iris[,5])
irisTargets <- decodeClassLabels(iris[,5], valTrue=0.9, valFalse=0.1)

iris <- splitForTrainingAndTest(irisValues, irisTargets, ratio=0.15)

#normalize data
iris <- normTrainingAndTestSet(iris)


#parameterGrid <- expand.grid(c(3,5,9,15), c(0.00316, 0.0147, 0.1))
parameterGrid <- expand.grid(c(3,5,9,15), c(0.00316, 0.0147, 0.1))

colnames(parameterGrid) <- c("nHidden", "learnRate")
rownames(parameterGrid) <- paste("nnet-", apply(parameterGrid, 1, function(x) {paste(x,sep="", collapse="-")}), sep="")

models <- apply(parameterGrid, 1, function(p) {
      
      mlp(iris$inputsTrain, iris$targetsTrain, size=p[1], learnFunc="Std_Backpropagation",
          learnFuncParams=c(p[2], 0.1), maxit=200, inputsTest=iris$inputsTest, 
          targetsTest=iris$targetsTest)
    })

par(mfrow=c(4,3))

for(modInd in 1:length(models)) {
  plotIterativeError(models[[modInd]], main=names(models)[modInd])
}

trainErrors <- data.frame(lapply(models, function(mod) {
          error <- sqrt(sum((mod$fitted.values - iris$targetsTrain)^2))
          error
        }))

testErrors <- data.frame(lapply(models, function(mod) {
      pred <- predict(mod,iris$inputsTest)
      error <- sqrt(sum((pred - iris$targetsTest)^2))
      error
    }))

t(trainErrors)
t(testErrors)

trainErrors[which(min(trainErrors) == trainErrors)]
testErrors[which(min(testErrors) == testErrors)]

model <- models[[which(min(testErrors) == testErrors)]]

model

summary(model)

