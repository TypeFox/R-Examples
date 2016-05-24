simrelplot <- function(obj, ncomp=min(obj$p,obj$n,20), ask=FALSE, print.cov=FALSE){

  def.par <- par(no.readonly = TRUE)
  if(!ask){
    layout(matrix(c(1,1,2,3),2,2, byrow=TRUE))
  }
  par(mar=c(5.1, 4.1, 4.1, 4.1))
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  dev.hold()
  plot(obj$beta, type="h", lwd=1, xlab="Variable number", ylab=expression(beta),
       main="True regression coefficients")
  abline(h=0, col=2)
  dev.flush()
  
  covs <- obj$Sigma[2:obj$p,1]
  covs <- abs(covs)/max(abs(covs))
  dev.hold()
    plot(obj$lambda[1:ncomp], type="h", lwd=2, col=1, 
       main="Relevant components plot",
       xlab="Components", ylab="Eigenvalue", axes=F, ylim=c(0,1))
    points(1:ncomp, covs[1:ncomp], type="p", pch=20, cex=2, col=2)
    axis(1)
    axis(2,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
    axis(4,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
    mtext("Scaled covariance",side=4, line=3, cex=0.85)
    box()
  dev.flush()
  
  X <- scale(obj$X, center=TRUE, scale=FALSE)
  Y <- scale(obj$Y, center=TRUE, scale=FALSE)
  svdres <- svd(X)
  eigval <- (svdres$d^2)/(svdres$d^2)[1]
  Z <- X%*%svdres$v
  covs <- cov(Y, Z)
  covs.sc <- abs(covs)/max(abs(covs))
  dev.hold()
    plot(eigval[1:ncomp], type="h", lwd=2, col=1, 
       main="Estimated relevant components plot",
       xlab="Components", ylab="Eigenvalue", axes=F, ylim=c(0,1))
    points(1:ncomp, covs.sc[1:ncomp], type="p", pch=20, cex=2, col=2)
    axis(1)
    axis(2,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
    axis(4,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
    mtext("Scaled covariance",side=4, line=3, cex=0.85)
    box()
  dev.flush() 
  if(print.cov){
    covs <- covs[1,1:min(ncomp, obj$p)]
    cat("Absolute values of estimated covariances\n")
    names(covs) <- paste("Component", 1:min(ncomp,obj$p),sep="")
    print(abs(round(covs,3)))
  }
  par(def.par)
}