plot.genlasso <- function(x, type=c("primal", "dual", "both"),
                          numbers=FALSE, vlines=TRUE, xlab, ylab, ...) {
  type = type[[1]]
  if (!(type %in% c("primal", "dual", "both"))) {
    stop("Invalid type, must be \"primal\", \"dual\", or \"both\".")
  }
  
  if (type=="both") {
    par(mfrow=c(1,2))
  }

  # If the path is trivial, just make up the lambda values
  if (length(x$lambda)==1 && x$lambda==Inf) {
    lambda = 1
    hit = NA
  }
  # Otherwise grab the lambda values
  else {
    lambda = x$lambda
    hit = x$hit
  }

  beta = x$beta
  u = x$u
  p = nrow(beta)
  m = nrow(u)
  
  if (x$completepath) {
    lambda = c(lambda,0)
    hit = c(hit,NA)
    beta = cbind(beta,x$bls)
    u = cbind(u,rep(0,m))
  }

  xlim = c(0,lambda[1]*1.10)
  if (numbers) d = 0.05*(xlim[2]-xlim[1])
  else d = 0
  
  if (type != "dual") {
    if (missing(xlab) || type=="both") xlab = expression(lambda)
    if (missing(ylab) || type=="both") ylab = expression(paste("Coordinates of ",beta))

    cols = rainbow(p)
    plot(c(),c(),xlim=c(xlim[1]-d,xlim[2]),ylim=range(beta),
         xlab=xlab,ylab=ylab,...)
    if (p!=0) {
      matplot(c(xlim[2],lambda),rbind(beta[,1],t(beta)),
              col=cols,type="l",lty=1,add=TRUE)
    }
    if (vlines) abline(v=lambda[hit],lty=2)
    if (vlines) abline(v=lambda[!hit],lty=2,col="red")
    if (p!=0 && numbers) text(xlim[1]-d,beta[,ncol(beta)],Seq(1,p))
  }

  if (type != "primal") {
    if (missing(xlab) || type=="both") xlab = expression(lambda)
    if (missing(ylab) || type=="both") ylab = "Coordinates of u"
          
    cols = rainbow(m)
    plot(c(),c(),xlim=c(xlim[1],xlim[2]+d),ylim=c(-xlim[2],xlim[2]),
         xlab=xlab,ylab=ylab, ...)
    lines(c(0,xlim[2]),c(0,xlim[2]),lwd=2)
    lines(c(0,xlim[2]),c(0,-xlim[2]),lwd=2)
    if (m!=0) {
      matplot(c(xlim[2],lambda),rbind(u[,1],t(u)),
              col=cols,type="l",lty=1,add=TRUE)
    }
    if (vlines) abline(v=lambda[hit],lty=2)
    if (vlines) abline(v=lambda[!hit],lty=2,col="red")
    if (m!=0 && numbers) text(xlim[2]+d,u[,1],Seq(1,m))
  }
}
