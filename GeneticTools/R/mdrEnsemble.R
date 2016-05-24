mdrEnsemble <- function(mdr, data=NULL, new.status=NULL, fold=NULL){

  if(is.null(data)) data <- mdr$X
  if(is.null(new.status)) new.status <- mdr$status

  if(is.null(fold)) fold <- mdr$fold
  top <- mdr$top
  t <- mdr$t
  oldX <- mdr$X
  oldStatus <- mdr$status

  if(fold==1){
    colOI <- matrix(mdr$mdr$topOneIndex[top:1],ncol=1)
  } else if(fold==2){
    colOI1 <- mdr$mdr$topTwoIndex1[top:1]
    colOI2 <- mdr$mdr$topTwoIndex2[top:1]
    colOI <- cbind(colOI1,colOI2)
  } else if(fold==3){
    colOI1 <- mdr$mdr$topThreeIndex1[top:1]
    colOI2 <- mdr$mdr$topThreeIndex2[top:1]
    colOI3 <- mdr$mdr$topThreeIndex3[top:1]
    colOI <- cbind(colOI1,colOI2,colOI3)
  } else if(fold==4){
    colOI1 <- mdr$mdr$topFourIndex1[top:1]
    colOI2 <- mdr$mdr$topFourIndex2[top:1]
    colOI3 <- mdr$mdr$topFourIndex3[top:1]
    colOI4 <- mdr$mdr$topFourIndex4[top:1]
    colOI <- cbind(colOI1,colOI2,colOI3,colOI4)
  }
  colOI <- colOI - 1
  res <- .Call("mdrEnsemble",X=data, fold=fold, status=new.status, t=t, top=colOI, oldStatus=oldStatus, oldX=oldX)

  if(fold==1){
    classifier <- res$classLableOne
  } else if (fold==2){
    classifier <- res$classLableTwo
  } else if(fold==3){
    classifier <- res$classLableThree
  } else if(fold==4){
    classifier <- res$classLableFour
  }
  tableClassifier <- apply(classifier,1,mdrTable)
  result <- c()
  for(i in 1:ncol(tableClassifier)){
    temp <- tableClassifier[,i]
    vote <- temp[2] / temp[1]
    result[i] <- as.numeric(vote>1)
    #result[i] <- vote
  }
  output <- list(result=result,cv=list(evalOne=res$evalOne,evalTwo=res$evalTwo,evalThree=res$evalThree,evalFour=res$evalFour))
  output
}

mdrTable <- function(x){
  values0 <- sum(x==0)
  values1 <- sum(x==1)
  c(values0,values1)
}