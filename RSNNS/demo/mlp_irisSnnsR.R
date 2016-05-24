library(RSNNS)

basePath <- ("./")

data(iris)

set.seed(2)

#normalize data
inputs <- normalizeData(iris[,1:4], "norm")

#outputs <- decodeClassLabels(iris[,5])
outputs <- decodeClassLabels(iris[,5], valTrue=0.9, valFalse=0.1)

numHiddenUnits <- 10

snnsObject <- SnnsRObjectFactory()

snnsObject$setLearnFunc('Quickprop')
snnsObject$setUpdateFunc('Topological_Order')
snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')

snnsObject$createNet(c(ncol(inputs),numHiddenUnits,ncol(outputs)), TRUE)

patset <- snnsObject$createPatSet(inputs, outputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$initializeNet(-1)

snnsObject$saveNet(paste(basePath,"mlp_irisSnnsR_untrained.net",sep=""),"mlp_irisSnnsR_untrained.net")

parameters <- c(0.2, 0, 0, 0, 0)
maxit <- 100

error <- vector()
for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
  error[i] <- res[[2]]
}

plot(error, type="l")

predictions <- snnsObject$predictCurrPatSet("output", c(0))

confusionMatrix(outputs,predictions)

snnsObject$saveNet(paste(basePath,"mlp_irisSnnsR.net",sep=""),"mlp_irisSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"mlp_irisSnnsR.pat",sep=""), patset$set_no)

