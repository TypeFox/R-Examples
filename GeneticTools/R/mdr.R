# Version: 30-11-2012, Daniel Fischer

# The option NAasValues is not working yet! The idea is here to create then 4x4 matrices and
# threat them as values in order to see, if they provide any information...

mdr <- function(X,status,fold=2,t=NULL,cv=0,cvp=0.75,top=20,NAasValues=TRUE,fix=NULL){

    N <- nrow(X)
    cvRes <- list()
    
    if(is.character(fix)){
      fix <- which((colnames(X)==fix)==TRUE)
    }

    ifelse(is.null(fix), fix <- 0,fix <- fix - 1)
    if(is.null(t)) t <- table(status)[2]/table(status)[1]
    if(nrow(X)!=length(status)) stop ("Nrow(X) / length(Status) mismatch!\n")

    res <- mdr.C(X=X,fold=fold,status=status,t=t,cv=0,cvp=cvp,top=top,na=as.numeric(NAasValues),fix)  

    if(cv>0){
      indices <- 1:nrow(X)
      for(i in 1:cv)
      { cvSub <- list()
	trainSet <- sample(indices,floor(length(status)*cvp))
	testSet <- indices[-trainSet]
	## ACHTUNG!!! BERECHNE HIER T NEU!!! IST NOCH FALSCH!!!!
	trainModel <- mdr.C(X=X[trainSet,],fold=fold,status=status[trainSet],t=t,cv=0,cvp=cvp,top=top,na=as.numeric(NAasValues),fix)
        tempModel <- list(mdr=trainModel,fold=fold,t=t,cv=0,top=top,fix=fix,X=X,status=status,cv=cvRes)
	for(foldRun in 1:fold)
	{
	  trainSet <- mdrEnsemble(tempModel,data=X[trainSet,],new.status=status[trainSet],fold=foldRun)$cv[[foldRun]]
          testSet <- mdrEnsemble(tempModel,data=X[testSet,],new.status=status[testSet],fold=foldRun)$cv[[foldRun]]
	  cvSub[[foldRun]] <- list(train=trainSet,test=testSet)
	}
        cvRes[[i]] <- cvSub
      }
    }

    result <- list(mdr=res,fold=fold,t=t,cv=cv,top=top,fix=fix,X=X,status=status,cvRes=cvRes)

    class(result) <- "mdr"
    result
} 

