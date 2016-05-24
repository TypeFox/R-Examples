library(RSNNS)

basePath <- ("./")

data(snnsData)
dataset <- snnsData$dlvq_ziff_100.pat

inputs <- dataset[,inputColumns(dataset)]
outputs <- dataset[,outputColumns(dataset)]

snnsObject <- SnnsRObjectFactory()

snnsObject$setLearnFunc('Dynamic_LVQ')
snnsObject$setUpdateFunc('Dynamic_LVQ')
snnsObject$setUnitDefaults(0,0,3,0,1,'Act_Logistic','Out_Identity')

snnsObject$createNet(c(ncol(inputs), 1), fullyConnectedFeedForward = FALSE)

patset <- snnsObject$createPatSet(inputs, outputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$initializeNet(c(1.0, -1.0, 0, 0, 0), "DLVQ_Weights")
snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$saveNet(paste(basePath,"dlvq_ziffSnnsR_untrained.net",sep=""),"dlvq_ziffSnnsR_untrained")

parameters <- c(0.03, 0.03, 10.0, 10.0, 0.0)

snnsObject$learnAllPatterns(parameters)

snnsObject$saveNet(paste(basePath,"dlvq_ziffSnnsR.net",sep=""),"dlvq_ziffSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"dlvq_ziffSnnsR.pat",sep=""), patset$set_no);

predictions <- snnsObject$predictCurrPatSet("output", parameters)

confusionMatrix(outputs, predictions)

mean(predictions - outputs)
