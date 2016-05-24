`new.penalized.pls` <-
function(ppls,Xtest,ytest=NULL){

  ncomp=ncol(ppls$coefficients)
  
  #if (!is.matrix(Xtest)) Xtest=matrix(Xtest,ncol=1)

  ypred<-Xtest%*%ppls$coefficients +rep(1,nrow(Xtest))%*%t(ppls$intercept)

  mse=NULL

  if (!is.null(ytest)){

    yres=ypred-ytest

    mse=apply(yres^2,2,mean)

  }

  return(list(ypred=ypred,mse=mse))

}
