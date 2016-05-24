plot.bayesm.mat=function(x,names, ...){
#
#  S3 method to print matrices of draws the object X is of class "bayesm.mat"
#  Modified after P. Rossi 2/07
#
  X=x
  if(mode(X) == "list") stop("list entered \n Possible Fixup: extract from list \n")
  if(mode(X) !="numeric") stop("Requires numeric argument \n")
  if(is.null(attributes(X)$dim)) X=as.matrix(X)
  nx=ncol(X)
  if(nx==1) par(mfrow=c(1,1)) 
  if(nx==2) par(mfrow=c(2,1))
  if(nx==3) par(mfrow=c(3,1))
  if(nx==4) par(mfrow=c(2,2))
  if(nx>=5) par(mfrow=c(3,2))
  if(missing(names)) {names=as.character(1:nx)}
  if(nx==1) par(mfrow=c(1,2)) 
  if(nx==2) par(mfrow=c(2,2))
  if(nx>=3) par(mfrow=c(3,2))
  for(index in 1:nx){
    plot(as.vector(X[,index]),xlab="",ylab="",main=names[index],type="l", col="blue")
    if(var(X[,index])>1.0e-20) {acf(as.vector(X[,index]),xlab="",ylab="",main="")}
    else { 
      plot.default(X[,index],xlab="",ylab="",type="n",main="No ACF Produced")}
    par(ask=TRUE) 
  }
  par(mfrow=c(1,1), ask=FALSE) 
  invisible()
}

