library(RSNNS)

data(snnsData)
trainData <- snnsData$artmap_train.pat
testData <- snnsData$artmap_test.pat

model <- artmap(trainData, nInputsTrain=70, nInputsTargets=5, nUnitsRecLayerTrain=50, nUnitsRecLayerTargets=26)
#updateFunc="ARTMAP_Synchronous"

model$fitted.values

predict(model, testData)

#summary(model)
