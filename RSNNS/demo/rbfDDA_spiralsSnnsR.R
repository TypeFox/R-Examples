library(RSNNS)

basePath <- ("./")

data(snnsData)
inputs <- snnsData$spirals.pat[,inputColumns(snnsData$spirals.pat)]
outputs <- snnsData$spirals.pat[,outputColumns(snnsData$spirals.pat)]

snnsObject <- SnnsRObjectFactory()

snnsObject$setLearnFunc('RBF-DDA')
snnsObject$setUpdateFunc('Topological_Order')
snnsObject$setUnitDefaults(0,0,1,0,1,'Act_Logistic','Out_Identity')

#snnsObject$setLearnFunc('RadialBasisLearning')
#snnsObject$setUpdateFunc('Topological_Order')
#snnsObject$setUnitDefaults(0,0,1,0,1,'Act_RBF_Gaussian','Out_Identity')
#snnsObject$setInitialisationFunc('RBF_Weights')
#snnsObject$createNet(c(2,45,2), TRUE)

snnsObject$createNet(c(2,2), fullyConnectedFeedForward = FALSE)

patset <- snnsObject$createPatSet(inputs, outputs)
snnsObject$setCurrPatSet(patset$set_no)

#snnsObject$initializeNet(c(1.0,  -1.0,  0.3,  1.0,  0.5) )
#snnsObject$initializeNet(0)
snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$saveNet(paste(basePath,"rbfDDA_spiralsSnnsR_untrained.net",sep=""),"rbfDDA_spiralsSnnsR_untrained")

parameters <- c(0.4, 0.2, 5)

res <- snnsObject$learnAllPatterns(parameters)

#maxit <- 100
#for(i in 1:maxit) {
#}

predictions <- snnsObject$predictCurrPatSet("output", c(0))

p <- encodeClassLabels(predictions, method="WTA", l=0, h=0)
t <- encodeClassLabels(outputs)

confusionMatrix(t,p)

snnsObject$saveNet(paste(basePath,"rbfDDA_spiralsSnnsR.net",sep=""),"rbfDDA_spiralsSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"rbfDDA_spiralsSnnsR.pat",sep=""), patset$set_no);









