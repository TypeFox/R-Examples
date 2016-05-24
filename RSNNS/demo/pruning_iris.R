library(RSNNS)

set.seed(2)

data(iris)

#shuffle the vector
df <- iris[sample(nrow(iris)),]


dfValues <- df[,1:4]
#dfTargets <- decodeClassLabels(df[,5])
dfTargets <- decodeClassLabels(df[,5], valTrue=0.9, valFalse=0.1)

df <- splitForTrainingAndTest(dfValues, dfTargets, ratio=0.15)

#normalize data
df <- normTrainingAndTestSet(df)

#df <- list()
#
##shuffle the vector
#iris <- iris[sample(nrow(iris)),]
#
##normalize data
#df$inputsTrain <- normalizeData(iris[,1:4], "norm")
#
##outputs <- decodeClassLabels(iris[,5])
#df$targetsTrain <- decodeClassLabels(iris[,5], valTrue=0.9, valFalse=0.1)




pruneFuncParams <- list(max_pr_error_increase = 10.0, pr_accepted_error = 1.0, no_of_pr_retrain_cycles = 1000, min_error_to_stop = 0.01,
                            init_matrix_value = 1e-6, input_pruning = TRUE, hidden_pruning = TRUE)
         
model <- mlp(df$inputsTrain, df$targetsTrain, size=c(10,10), learnFunc="Std_Backpropagation", learnFuncParams=c(0.2, 0, 0, 0, 0), 
             pruneFunc="OptimalBrainSurgeon", pruneFuncParams=pruneFuncParams, maxit=100000, inputsTest=df$inputsTest, 
             targetsTest=df$targetsTest)

par(mfrow=c(2,2))

plotIterativeError(model)

predictions <- predict(model,df$inputsTest)

plotRegressionError(df$targetsTest[,2], predictions[,2])

confusionMatrix(df$targetsTrain,fitted.values(model))
confusionMatrix(df$targetsTest,predictions)

plotROC(fitted.values(model)[,2], df$targetsTrain[,2])
plotROC(predictions[,2], df$targetsTest[,2])

#confusion matrix with 402040-method
confusionMatrix(df$targetsTrain, encodeClassLabels(fitted.values(model),method="402040", l=0.4, h=0.6))

model

weightMatrix(model)

summary(model)

extractNetInfo(model)
