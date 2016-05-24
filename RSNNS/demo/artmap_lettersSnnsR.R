
library(RSNNS)

basePath <- ("./")
#basePath <- ("/home/bergmeir/")

data(snnsData)
trainData <- snnsData$artmap_train.pat
testData <- snnsData$artmap_test.pat

inputs <- trainData

snnsObject <- SnnsRObjectFactory()

#bn_artmap_createNet(int f1aUnits, int f1aRows, int f2aUnits, 
#    int f2aRows, int f1bUnits, int f1bRows, 
#    int f2bUnits, int f2bRows)

snnsObject$setInitialisationFunc('ARTMAP_Weights')
snnsObject$setLearnFunc('ARTMAP')
snnsObject$setUpdateFunc('ARTMAP_Stable')
#snnsObject$setUpdateFunc('ARTMAP_Synchronous')

#snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')

snnsObject$artmap_createNet(70,14,50,14,5,1,26,6)

snnsObject$saveNet(paste(basePath,"artmap_SnnsR_untrained.net",sep=""),"artmap_SnnsR_untrained")

patset <- snnsObject$createPatSet(inputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$initializeNet(c(1.0, 1.0, 1.0, 1.0, 0.0))

snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

parameters <- c(0.8, 1.0, 1.0, 0, 0)

snnsObject$learnAllPatterns(parameters)

snnsObject$saveNet(paste(basePath,"artmap_SnnsR.net",sep=""),"artmap_SnnsR")
snnsObject$saveNewPatterns(paste(basePath,"artmap_SnnsR.pat",sep=""), patset$set_no);

patsetTest <- snnsObject$createPatSet(testData)
snnsObject$setCurrPatSet(patsetTest$set_no)

outputs <- snnsObject$predictCurrPatSet("artmap", parameters)

outputs

