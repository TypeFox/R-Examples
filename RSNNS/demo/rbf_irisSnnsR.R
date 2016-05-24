library(RSNNS)

basePath <- ("./")

data(iris)

set.seed(2)

#normalize data
inputs <- normalizeData(iris[,1:4], "norm")

#outputs <- decodeClassLabels(iris[,5])
outputs <- decodeClassLabels(iris[,5], valTrue=0.9, valFalse=0.1)

numHiddenUnits <- 40

snnsObject <- SnnsRObjectFactory()

snnsObject$setLearnFunc('RadialBasisLearning')
snnsObject$setUpdateFunc('Topological_Order')
snnsObject$setUnitDefaults(0,0,1,0,1,'Act_RBF_Gaussian','Out_Identity')
snnsObject$createNet(c(ncol(inputs),numHiddenUnits,ncol(outputs)), TRUE)

snnsObject$setTTypeUnitsActFunc("UNIT_INPUT", "Act_Identity")
snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_RBF_Gaussian")
#snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", "Act_IdentityPlusBias")
snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", "Act_Logistic")

patset <- snnsObject$createPatSet(inputs, outputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$initializeNet(c(0,0,0), "RBF_Weights_Kohonen")
snnsObject$initializeNet(c(-4, 4, 0, 0.2, 0.05), "RBF_Weights")
#snnsObject$initializeNet(c(0, 1, 0, 0.02, 0.01), "RBF_Weights")
#snnsObject$learnAllPatterns(c(0,0,0,0.1,0.9))
#snnsObject$testAllPatterns(c(0,0,0,0.1,0.9))


snnsObject$saveNet(paste(basePath,"rbf_irisSnnsR_untrained.net",sep=""),"rbf_irisSnnsR_untrained.net")

#parameters <- c(0, 0, 0, 0.1, 0.9)
parameters <- c(0, 0, 0.01, 0.1, 0.8)
maxit <- 1000

error <- vector()
for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
  error[i] <- res[[2]]
}

plot(error, type="l")

predictions <- snnsObject$predictCurrPatSet("output", c(0))

confusionMatrix(outputs,predictions)

snnsObject$saveNet(paste(basePath,"rbf_irisSnnsR.net",sep=""),"rbf_irisSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"rbf_irisSnnsR.pat",sep=""), patset$set_no)

