#' Plot results of an AME object
#' 
#' A set of plots summarizing the MCMC routine for an AME fit, as well as some
#' posterior predictive checks.
#' 
#' 
#' @param x the result of fitting an AME model
#' @param ... additional parameters (not used)
#' @return a series of plots
#' @author Peter Hoff
#' @S3method plot ame
plot.ame <-
function(x, ...)
{  
  fit<-x 
  require(amen) 

  gof<-1*(nrow(fit$GOF)>1) 
  par(mfrow=c(1+2*gof,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

  mVC<-apply(fit$VC,2,median)
  matplot(fit$VC,type="l",lty=1,ylab="VC")
  abline(h=mVC,col=1:length(mVC) )

  if(ncol(fit$BETA)>0) 
  { 
    mBETA<-apply(fit$BETA,2,median)
    matplot(fit$BETA,type="l",lty=1,col=1:length(mBETA),ylab="BETA")
    abline(h=mBETA,col=1:length(mBETA) )
    abline(h=0,col="gray")
  }

  if(gof)
  {
    for(k in 1:4)
    {
      hist(fit$GOF[-1,k],xlim=range(fit$GOF[,k]),main="",prob=TRUE,
           xlab=colnames(fit$GOF)[k],col="lightblue",ylab="",yaxt="n")
           abline(v=fit$GOF[1,k],col="red")
    }
  }
}

