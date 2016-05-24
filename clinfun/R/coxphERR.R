coxphERR <- function(phfit, ngamma=NULL) {
  if (class(phfit) != "coxph") stop("phfit shoud be coxph class object")
  if (is.null(phfit$x)) stop("coxph should have been called with x=TRUE option")
  ss <- phfit$n
  covec <- phfit$x-matrix(rep(apply(phfit$x,2,mean),ss), nrow=ss, byrow=T)
  if (is.null(ngamma)) ngamma <- 1:ncol(covec)
  covec.gamma <- covec[,ngamma,drop=FALSE]
  covec.beta <- covec[,-ngamma,drop=FALSE]
  xf <- exp(covec%*%phfit$coef)
  xf.gamma <- exp(covec.gamma%*%phfit$coef[ngamma])
  ERR <- (log(mean(xf.gamma)))/(0.577215665+(log(mean(xf))))

  g <- log(mean(xf))
  h <- log(mean(xf.gamma))
  dgb <- (t(xf)%*%covec.beta)/sum(xf)
  dga <- (t(xf)%*%covec.gamma)/sum(xf)
  dha <- (t(xf.gamma)%*%covec.gamma)/sum(xf.gamma)
  dxfb <- (-h*dgb)/((0.577215665+g)^2)
  dxfa <- (((0.577215665+g)*dha)-(h*dga))/((0.577215665+g)^2)
  dxf <-cbind(dxfa,dxfb)
  newcovmat11 <- phfit$var[ngamma,ngamma,drop=FALSE]
  newcovmat12 <- phfit$var[ngamma,-ngamma,drop=FALSE]
  newcovmat21 <- phfit$var[-ngamma,ngamma,drop=FALSE]
  newcovmat22 <- phfit$var[-ngamma,-ngamma,drop=FALSE]
  newcovmat <- cbind(rbind(newcovmat11,newcovmat21),rbind(newcovmat12,newcovmat22))
  vr <- var(xf.gamma)/ss
  vs <- var(xf)/ss
  cvrs <- cov(xf.gamma,xf)/ss
  
  dr <- 1/((mean(xf.gamma))*(0.577215665+log(mean(xf))))
  ds <- (-log(mean(xf.gamma)))/(mean(xf)*((0.577215665+log(mean(xf)))^2))
  dd <- cbind(dr,ds)
  newmat <- cbind(rbind(vr,cvrs),rbind(cvrs,vs))
  se.ERR <- sqrt((dd%*%newmat%*%t(dd)) + (dxf%*%newcovmat%*%t(dxf)))

  out <- c(ERR, se.ERR)
  names(out) <- c("ERR","se.ERR")
  out
}
