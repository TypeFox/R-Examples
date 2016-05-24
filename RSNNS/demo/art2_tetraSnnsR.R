library(RSNNS)

basePath <- ("./")

data(snnsData)
dataset <- snnsData$art2_tetra_med.pat

inputs <- dataset

snnsObject <- SnnsRObjectFactory()

snnsObject$setUnitDefaults(0,0,3,0,1,'Act_Logistic','Out_Identity')

snnsObject$art2_createNet(3,3,5,5)

snnsObject$setInitialisationFunc('ART2_Weights')
snnsObject$setLearnFunc('ART2')
snnsObject$setUpdateFunc('ART2_Stable')

patset <- snnsObject$createPatSet(inputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$initializeNet(c(0.9, 2.0, 0, 0, 0))
snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$saveNet(paste(basePath,"art2_tetraSnnsR_untrained.net",sep=""),"art2_tetraSnnsR_untrained")

parameters <- c(0.99, 20.0, 20.0, 0.1, 0.0)
maxit <- 100

for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
  if(res[[1]] != 0) print(paste("An error occured at iteration ", i, " : ", res, sep=""))
}

snnsObject$saveNet(paste(basePath,"art2_tetraSnnsR.net",sep=""),"art2_tetraSnnsR")
snnsObject$saveNewPatterns(paste(basePath,"art2_tetraSnnsR_train.pat",sep=""), patset$set_no);

#testdata <- snnsData$art2_tetra_med.pat
#testPatset <- snnsObject$createPatSet(testdata)
#snnsObject$setCurrPatSet(testPatset$set_no)
#snnsObject$saveNewPatterns(paste(basePath,"art2_tetraSnnsR_test.pat",sep=""), patset$set_no);
outputs <- snnsObject$predictCurrPatSet("art2", parameters)
outputs

encodeClassLabels(outputs)

library(scatterplot3d)
scatterplot3d(inputs, pch=encodeClassLabels(outputs))

##------------------------------------------------------------------------------------------------
## using the artui interface of SNNS instead of the normal functions
## however, the artui interface is currently not included in the wrapping as it seems to be rather
## unstable and the normal interface can be used instead..
#noOfPatterns <- snnsObject$getNoOfPatterns()
#updateFuncParams <- parameters
#
#print(snnsObject$artui_getN())
#print(snnsObject$artui_getM())
#
#classNumbers <- NULL
#for(currentPattern in 1:noOfPatterns) {
#  snnsObject$setPatternNo(currentPattern)
#  snnsObject$showPattern(resolveSnnsRDefine("patternUpdateModes","OUTPUT_NOTHING"))
#  snnsObject$updateNet(updateFuncParams)
#  
#  print(snnsObject$artui_getClassifiedStatus()$status)
#  classNumbers <- c(classNumbers, snnsObject$artui_getClassNo()$class_no)
#  
#}
#classNumbers
