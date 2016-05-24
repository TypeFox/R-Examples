library(RSNNS)

basePath <- "./"

data(snnsData)
inputs <- snnsData$art1_letters.pat

snnsObject <- SnnsRObjectFactory()

snnsObject$setInitialisationFunc('ART1_Weights')
snnsObject$setLearnFunc('ART1')
snnsObject$setUpdateFunc('ART1_Synchronous')
snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')

snnsObject$art1_createNet(35,7,26,7)

patset <- snnsObject$createPatSet(inputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$initializeNet(c(1.0, 1.0))
snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$saveNet(paste(basePath,"art1_lettersSnnsR_untrained.net",sep=""),"art1_lettersSnnsR_untrained")

parameters <- c(0.9, 0, 0)
maxit <- 100

for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
  if(res[[1]] != 0) print(paste("An error occured at iteration ", i, " : ", res, sep=""))
}

snnsObject$saveNet(paste(basePath,"art1_lettersSnnsR.net",sep=""),"art1_lettersSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"art1_lettersSnnsR.pat",sep=""), patset$set_no);

outputs <- snnsObject$predictCurrPatSet("art1", parameters)
encodeClassLabels(outputs)

inputMaps <- matrixToActMapList(inputs, nrow=7)
par(mfrow=c(3,3))
for (i in 1:9) plotActMap(inputMaps[[i]])
