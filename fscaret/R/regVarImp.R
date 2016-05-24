regVarImp <- function(model, xTrain, yTrain, xTest, fitControl, myTimeLimit, no.cores, lk_col, supress.output, mySystem){

resultVarImpListCombREG <- NULL
resultVarImpListCombREG <- list()

regVarImpFuncRES <- NULL
regVarImpFuncRES <- list()

myTimeLimitSet <- myTimeLimit
fitControlSet <- fitControl
lk_col <- lk_col
supress.output <- supress.output
mySystem <- mySystem
no.cores <- no.cores


regVarPred <- function(funcRegPred) {

#Print out function names
cat("----------------------------------------\n")
cat("Calculating: ", funcRegPred,"\n")
cat("----------------------------------------\n")


outfile<-paste("./",(lk_col-1),"in_default_REGControl_", paste(funcRegPred),".RData",sep="")

outfileImp<-paste("./",(lk_col-1),"in_default_REGControl_VarImp_", paste(funcRegPred),".txt",sep="")


#start feature selection method
timer1 <- proc.time()

if(mySystem!="windows"){
  if(supress.output==TRUE){

    # Supress output
    sink("/dev/null")
      res <- invisible(try(timeout(train(xTrain,yTrain, method=funcRegPred, trControl=fitControlSet),seconds=myTimeLimitSet),silent=TRUE))
    sink()

    } else {

      res <- invisible(try(timeout(train(xTrain,yTrain, method=funcRegPred, trControl=fitControlSet),seconds=myTimeLimitSet),silent=TRUE))
    }
    
} else {

    if(supress.output==TRUE){

    # Supress output
    sink("NUL")
    res <- try(train(xTrain,yTrain, method=funcRegPred, trControl=fitControlSet),silent=TRUE)
    sink()

    } else {

    res <- try(train(xTrain,yTrain, method=funcRegPred, trControl=fitControlSet),silent=TRUE)

    }

}

timer2 <- proc.time() - timer1 

variableImportanceRes <- try(varImp(res$finalModel),silent=TRUE)

resultVarImpListCombREG[funcRegPred] <- try(list(variableImportanceRes),silent=TRUE)

      cat("----------------------------------------\n")
      cat("",funcRegPred,"\n")
      cat("----------------------------------------\n")
      cat("Elapsed time: ",timer2,"\n")
#       cat("Variable importance: \n")

if((class(res) != "try-error")&&(class(variableImportanceRes) != "try-error")){

# save results
try(save(res, file=outfile),silent=TRUE)

try(write.table(variableImportanceRes, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp),silent=TRUE)

  } else if(class(res) != "try-error"){

      yPred <- try(predict(res,xTest),silent=TRUE)

    if(class(yPred)!="try-error"){

      variableImportanceRes <- try(filterVarImp(xTest,yPred,nonpara=TRUE),silent=TRUE)
      
	if(class(variableImportanceRes)!="try-error") {

	  resultVarImpListCombREG[funcRegPred] <- try(list(variableImportanceRes),silent=TRUE)
	  try(write.table(variableImportanceRes,col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp))
	  
	} else if(class(variableImportanceRes)=="try-error") {
	
	  print("Predicting variable importance (first try) has failed!")
	  resultVarImpListCombREG[funcRegPred] <- try(list(NA),silent=TRUE)
	
	}

      }
      
	if((class(res)!="try-error") && (class(variableImportanceRes) != "try-error")){
	  # save results
	    try(save(res, file=outfile),silent=TRUE)

	    try(write.table(variableImportanceRes,col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp))
      
	}

    if ((class(res)!="try-error") && (class(variableImportanceRes) == "try-error")){
    
    variableImportanceRes <- try(varImp(res),silent=TRUE)

    resultVarImpListCombREG[funcRegPred] <- try(list(variableImportanceRes$importance),silent=TRUE)
    
      # save results
      try(save(res, file=outfile),silent=TRUE)

      try(write.table(variableImportanceRes$importance,col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp))

    }

  } else if (class(variableImportanceRes) == "try-error"){

      variableImportanceRes <- try(varImp(res),silent=TRUE)

      if(class(variableImportanceRes)!="try-error") {

	resultVarImpListCombREG[funcRegPred] <- try(list(variableImportanceRes$importance),silent=TRUE)

	# save results
	try(save(res, file=outfile),silent=TRUE)

	try(write.table(variableImportanceRes$importance,col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp))

      } else if(class(variableImportanceRes)!="try-error"){

	    yPred <- try(predict(res,xTest),silent=TRUE)

	  if(class(yPred)!="try-error"){

	    variableImportanceRes <- try(filterVarImp(xTest,yPred,nonpara=TRUE),silent=TRUE)

	      if(class(variableImportanceRes)!="try-error") {

		resultVarImpListCombREG[funcRegPred] <- try(list(variableImportanceRes),silent=TRUE)

		try(write.table(variableImportanceRes,col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp))

	      } else if(class(variableImportanceRes)=="try-error") {

	      print("Predicting variable importance (second try) has failed!")
	      resultVarImpListCombREG[funcRegPred] <- try(list(NA),silent=TRUE)

	  }

	}
      
      } else if(class(res)=="try-error"){
      
	print("Building model has failed or reached timelimit!")
	resultVarImpListCombREG[funcRegPred] <- try(list(NA),silent=TRUE)
      }
      
  } else {
  
  print("Predicting variable importance (third try) has failed!")
  resultVarImpListCombREG[funcRegPred] <- try(list(NA),silent=TRUE)
  
  }

  # Last check to get all possible varImp  
    tmpSum <- sum(variableImportanceRes[,1])

  if(tmpSum == 0){
  
    variableImportanceRes <- try(varImp(res),silent=TRUE)
  
    resultVarImpListCombREG[funcRegPred] <- try(list(variableImportanceRes$importance),silent=TRUE)
  
    try(write.table(variableImportanceRes$importance,col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t", file=outfileImp))        
  }

  
  if(class(variableImportanceRes)!="try-error"){

  # Print out variable importance
  #   try(print(variableImportanceRes),silent=TRUE)
  cat("Variable importance has been calculated! \n")
  
  } else if(class(variableImportanceRes)=="try-error"){
  
  print("Predicting variable importance (fourth try) has failed!")
  resultVarImpListCombREG[funcRegPred] <- try(list(NA),silent=TRUE)
  
  }

  if((file.exists(paste(outfile))) || file.exists(outfileImp)){
#   general check for files RData and VarImp.txt
    if(((file.exists(paste(outfile)))==FALSE) || (file.exists(outfileImp)==FALSE)){
      try(file.remove(paste(outfile)))
      try(file.remove(paste(outfileImp)))
    } else {
    
    print("RData and VarImp.txt files exists!")
    cat("Variable importance: \n")
    try(print(variableImportanceRes),silent=TRUE)
    
    }
    
  }
  
}


# resultVarImpListCombREG[model] <- mclapply(model,regVarPred, mc.preschedule=TRUE, mc.cores=no.cores, mc.set.seed=TRUE, mc.cleanup=TRUE, mc.allow.recursive=TRUE)

if (Sys.info()[1] == "Windows"){

# Windows parallel implementation

# Spawn child processes using fork()
cl <- makeCluster(no.cores)

# Export objects to the cluster
clusterExport(
              cl=cl, 
              varlist=c("myTimeLimitSet", "fitControlSet", "lk_col", "supress.output",
                        "mySystem", "no.cores", "xTrain","yTrain", "funcRegPred", "fitControlSet")
              ,envir=environment())

# Run function
resultVarImpListCombREG[model] <- parLapply(cl, model, regVarPred)

# Stop cluster and kill child processes
stopCluster(cl)

} else {
  
# POSIX OS parallel implementation  
  resultVarImpListCombREG[model] <- mclapply(model,regVarPred,
                                             mc.preschedule=FALSE, mc.cores=no.cores, mc.cleanup=FALSE)  
  
  
}


# Return variable importance or NULL

for(i in 1:length(resultVarImpListCombREG)){
  
  if(class(resultVarImpListCombREG[i])=="try-error"){
    
    resultVarImpListCombREG[i] <- NULL
    
  }
  
}

return(resultVarImpListCombREG)

}