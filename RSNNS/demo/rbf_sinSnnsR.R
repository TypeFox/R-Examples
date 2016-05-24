library(RSNNS)

basePath <- ("./")

#inputs <- as.matrix(seq(0,1,0.01))
#outputs <- as.matrix(sin(x))

set.seed(2)

inputs <- as.matrix(seq(0,10,0.1))
outputs <- as.matrix(sin(inputs) + runif(inputs*0.2))
#outputs <- as.matrix(sin(inputs))
outputs <- normalizeData(outputs, "0_1")

numHiddenUnits <- 40

snnsObject <- SnnsRObjectFactory()

snnsObject$setLearnFunc('RadialBasisLearning')
snnsObject$setUpdateFunc('Topological_Order')
snnsObject$setUnitDefaults(0,0,1,0,1,'Act_RBF_Gaussian','Out_Identity')
snnsObject$createNet(c(ncol(inputs),numHiddenUnits,ncol(outputs)), TRUE)

snnsObject$setTTypeUnitsActFunc("UNIT_INPUT", "Act_Identity")
snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_RBF_Gaussian")
snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", "Act_IdentityPlusBias")
#snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", "Act_Logistic")

patset <- snnsObject$createPatSet(inputs, outputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$initializeNet(c(0,0,0,0,0), "RBF_Weights_Kohonen")

snnsObject$initializeNet(c(0, 1, 0, 0.01, 0.01), "RBF_Weights")
#snnsObject$initializeNet(c(-4, 4, 0, 0.01, 0.01), "RBF_Weights")

snnsObject$saveNet(paste(basePath,"rbf_sinSnnsR_untrained.net",sep=""),"rbf_sinSnnsR_untrained.net")

parameters <- c(1e-8, 0, 1e-8, 0.1, 0.8)
maxit <- 1000

error <- vector()
for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
  error[i] <- res[[2]]
}

par(mfrow=c(2,1))

plot(error, type="l")

predictions <- snnsObject$predictCurrPatSet("output", c(0))

plot(inputs, outputs)
lines(inputs, predictions, col="green")

snnsObject$saveNet(paste(basePath,"rbf_sinSnnsR.net",sep=""),"rbf_sinSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"rbf_sinSnnsR.pat",sep=""), patset$set_no);

