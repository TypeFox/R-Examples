basePath <- ("./")

library(RSNNS)

mySnnsObject <- SnnsRObjectFactory()

mySnnsObject$setLearnFunc('Quickprop')
mySnnsObject$setUpdateFunc('Topological_Order')
mySnnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')

print("Creating input layer")

inputs <- vector()
for(i in 1:9) {
 
  num <- mySnnsObject$createDefaultUnit()
  inputs[i] <- num
  
  mySnnsObject$setUnitName(num,paste("Input_",i,sep=""))
  mySnnsObject$setUnitPosition(num, i, 0, 0)
}

print("Creating hidden layer")

hidden <- vector()
for(i in 1:4) {
  
  num <- mySnnsObject$createDefaultUnit()
  hidden[i] <- num
  
  mySnnsObject$setUnitName(num,paste("Hidden_",i,sep=""))
  
  #HIDDEN
  mySnnsObject$setUnitTType(num,3)
  
  mySnnsObject$setUnitPosition(num, i+3, 2, 0)
  mySnnsObject$setCurrentUnit(num)
  for(j in inputs)  {
    mySnnsObject$createLink(j,0);
  };  
}

print("Creating output layer")

outputs <- vector()
for(i in 1:9) {
 
  num <- mySnnsObject$createDefaultUnit()
  outputs[i] <- num
  
  mySnnsObject$setUnitName(num,paste("Output_",i,sep=""))
  
  #OUTPUT
  mySnnsObject$setUnitTType(num,2)
  
  mySnnsObject$setUnitPosition(num, i, 4, 0)
  mySnnsObject$setCurrentUnit(num)
  for(j in hidden)  {
    mySnnsObject$createLink(j,0);
  };  
}

print("Creating patterns")

mySnnsObject$deleteAllPatterns()

patset <- mySnnsObject$allocNewPatternSet()

for(unum in 1:9)  {
  for(curnum in 1:9)  {
    mySnnsObject$setUnitActivation(inputs[curnum], as.numeric(curnum == unum));
    mySnnsObject$setUnitActivation(outputs[curnum], as.numeric(curnum == unum));
  }
  mySnnsObject$newPattern();
}

mySnnsObject$initializeNet(-1)
mySnnsObject$shufflePatterns(TRUE)
mySnnsObject$DefTrainSubPat()

mySnnsObject$saveNet(paste(basePath,"encoderSnnsCLib_untrained.net",sep=""),"encoderSnnsCLib_untrained")

error <- vector()
for(i in 1:500) {
  res <- mySnnsObject$learnAllPatterns(c(0.2, 0, 0, 0, 0))
  if(res[[1]] != 0) print(paste("An error occured at iteration ", i, " : ", res, sep=""))
  error[i] <- res[[2]]
}

error
plot(error, type="l")

mySnnsObject$getCompleteWeightMatrix()

mySnnsObject$saveNet(paste(basePath,"encoderSnnsCLib.net",sep=""),"encoderSnnsCLib")
mySnnsObject$saveNewPatterns(paste(basePath,"encoderSnnsCLib.pat",sep=""), patset$set_no);
