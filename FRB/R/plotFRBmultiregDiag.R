
diagplot <- function(x,...) UseMethod("diagplot")

diagplot.FRBmultireg <- function(x, Xdist = TRUE, ...) {

FRBres <- x

Y <- FRBres$Y
X <- FRBres$X
n <- nrow(X)
q <- ncol(Y)
# compute the residual distances:
resids <- Y - X %*% FRBres$coefficients
#if (!is.null(FRBres$int)) {   # for GS-estimates
#    resids <- resids - as.matrix(rep(1,n))%*%FRBres$int 
#    X <- cbind(rep(1,n), X)
#}
residsD <- sqrt(mahalanobis(resids, rep(0,q), FRBres$Sigma))

if (Xdist) {
  # now perform multivariate location/scatter estimation via the same estimator
  interceptdetection <- apply(X==1, 2, all)
  withoutind <- (1:ncol(X))[interceptdetection==FALSE]
  XwI <- as.matrix(X[,withoutind])

  if (FRBres$method$est=="S")
    ests <- Sest_multireg(as.matrix(rep(1,n)), XwI, bdp=FRBres$method$bdp, control=FRBres$control)
  else if (FRBres$method$est=="GS")
    ests <- GSest_multireg(as.matrix(rep(1,n)), XwI, bdp=FRBres$method$bdp, control=FRBres$control)
  else
    ests <- MMest_multireg(as.matrix(rep(1,n)), XwI, control=FRBres$control)

  Xloc <- ests$coefficients
  Xscatter <- ests$Sigma
  leverageD <- sqrt(mahalanobis(XwI, Xloc, Xscatter))

  xlimm <- range(leverageD)
  xlimm[2] <- xlimm[2] + (xlimm[2]-xlimm[1])*0.03
  plot(leverageD, residsD, xlab="Distance in X-space", ylab="Residual distance", cex.lab=1.3, pch=20, xlim=xlimm)
  title("Diagnostic plot based on robust estimates")
  abline(h=sqrt(qchisq(.975, q)))
  abline(v=sqrt(qchisq(.975, ncol(XwI))))

  # identify outliers
  inds <- 1:n
  badinds <- inds[residsD > sqrt(qchisq(.975, q))]
  if (length(badinds) <= 20) {
    for (j in badinds)
      text(leverageD[j], residsD[j], j, pos=4, cex=0.8, offset=0.3)
  }
}
else {
  plot(1:n, residsD, xlab="Index", ylab="Residual distance", cex.lab=1.3, pch=20)
  title("Diagnostic plot based on robust estimates")
  abline(h=sqrt(qchisq(.975, q)))

  # identify outliers
  inds <- 1:n
  badinds <- inds[residsD > sqrt(qchisq(.975, q))]
  if (length(badinds) <= 20) {
    for (j in badinds)
      text(inds[j], residsD[j], j, pos=4, cex=0.8, offset=0.3)
  }
}

}

