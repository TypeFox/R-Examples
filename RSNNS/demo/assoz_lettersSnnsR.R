library(RSNNS)

basePath <- ("./")

data(snnsData)

patterns <- snnsData$art1_letters.pat

snnsObject <- SnnsRObjectFactory()

snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Identity','Out_Identity')

# this function also sets learning, initialization and update function..
snnsObject$assoz_createNet(5,7)

patset <- snnsObject$createPatSet(patterns)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$initializeNet(c(1.0, -1.0))
snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$saveNet(paste(basePath,"assoz_lettersSnnsR_untrained.net",sep=""),"assoz_lettersSnnsR_untrained")

parameters <- c(0.01, 100, 0.0, 0.0, 0.0)
maxit <- 100

for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
}

snnsObject$saveNet(paste(basePath,"assoz_lettersSnnsR.net",sep=""),"assoz_lettersSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"assoz_lettersSnnsR.pat",sep=""), patset$set_no);

outputs <- snnsObject$predictCurrPatSet("assoz", c(50,0,0,0,0))
outputMaps <- matrixToActMapList(outputs, nrow=7)

par(mfrow=c(3,3))
for (i in 1:9) plotActMap(outputMaps[[i]])
