plot.bayesm.mat=function(x,names,burnin=trunc(.1*nrow(X)),tvalues,TRACEPLOT=TRUE,DEN=TRUE,INT=TRUE,
      CHECK_NDRAWS=TRUE,...){
#
#  S3 method to print matrices of draws the object X is of class "bayesm.mat"
#
#     P. Rossi 2/07
#
  X=x
  if(mode(X) == "list") stop("list entered \n Possible Fixup: extract from list \n")
  if(mode(X) !="numeric") stop("Requires numeric argument \n")
  op=par(no.readonly=TRUE)
  on.exit(par(op))
  on.exit(devAskNewPage(FALSE),add=TRUE)
  if(is.null(attributes(X)$dim)) X=as.matrix(X)
  nx=ncol(X)
  if(nrow(X) < 100 & CHECK_NDRAWS) {cat("fewer than 100 draws submitted \n"); return(invisible())}
  if(!missing(tvalues)){ 
        if(mode(tvalues) !="numeric") {stop("tvalues must be a numeric vector \n")} 
      else 
        {if(length(tvalues)!=nx) stop("tvalues are wrong length \n")}
  }
  if(nx==1) par(mfrow=c(1,1)) 
  if(nx==2) par(mfrow=c(2,1))
  if(nx==3) par(mfrow=c(3,1))
  if(nx==4) par(mfrow=c(2,2))
  if(nx>=5) par(mfrow=c(3,2))

  if(missing(names)) {names=as.character(1:nx)}
  if (DEN) ylabtxt="density" else ylabtxt="freq"
  devAskNewPage(TRUE)
  for(index in 1:nx){
     hist(X[(burnin+1):nrow(X),index],xlab="",ylab=ylabtxt,main=names[index],freq=!DEN,col="magenta",...)
     if(!missing(tvalues)) abline(v=tvalues[index],lwd=2,col="blue")
     if(INT){
     quants=quantile(X[(burnin+1):nrow(X),index],prob=c(.025,.975))
     mean=mean(X[(burnin+1):nrow(X),index])
     semean=numEff(X[(burnin+1):nrow(X),index])$stderr
     text(quants[1],0,"|",cex=3.0,col="green")
     text(quants[2],0,"|",cex=3.0,col="green")
     text(mean,0,"|",cex=3.0,col="red")
     text(mean-2*semean,0,"|",cex=2,col="yellow")
     text(mean+2*semean,0,"|",cex=2,col="yellow")
     }
  }
  if(TRACEPLOT){
     if(nx==1) par(mfrow=c(1,2)) 
     if(nx==2) par(mfrow=c(2,2))
     if(nx>=3) par(mfrow=c(3,2))
     for(index in 1:nx){
        plot(as.vector(X[,index]),xlab="",ylab="",main=names[index],type="l",col="red")
        if(!missing(tvalues)) abline(h=tvalues[index],lwd=2,col="blue")
        if(var(X[,index])>1.0e-20) {acf(as.vector(X[,index]),xlab="",ylab="",main="")}
        else 
            {plot.default(X[,index],xlab="",ylab="",type="n",main="No ACF Produced")}
     }
   }
  invisible()
}

