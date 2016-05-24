##  
##  coded by E. S. Venkatraman based on the R code by Gonen 02/08/2007
##

coxphCPE <- function(phfit) {
  if (class(phfit) != "coxph") stop("phfit shoud be coxph class object")
  n <- phfit$n
  betahat <- phfit$coefficients
  p <- length(phfit$coefficients)
  vbetahat <- phfit$var
  xbeta <- phfit$linear.predictors
  if (getRversion() <= "2.9.1") {
    xMat <- as.matrix(model.matrix(phfit)[,-1])
  } else {
    xMat <- as.matrix(model.matrix(phfit))
  }
  bw <- 0.5*sd(xbeta)*(n^(-1/3))
  zzz <- .Fortran("cpesub",
                  as.integer(n),
                  as.integer(p),
                  as.double(xMat),
                  as.double(xbeta),
                  as.double(bw),
                  CPE=double(1),
                  CPEsmooth=double(1),
                  varDeriv=double(p),
                  uRowSum=double(n),
                  uSSQ=double(1),
                  PACKAGE="clinfun")
  CPE <- 2*zzz$CPE/(n*(n-1))
  CPEsmooth <- 2*zzz$CPEsmooth/(n*(n-1))
  varTerm1 <-  4*(sum((zzz$uRowSum+rep(0.5,n)-n*CPEsmooth)^2) - (zzz$uSSQ + n/4 - n*CPEsmooth - n*(n-2)*CPEsmooth^2))/(n*(n-1))^2
  varDeriv <- 2*zzz$varDeriv/(n*(n-1))
  varTerm2 <- t(varDeriv)%*%vbetahat%*%varDeriv
  varCPE <- varTerm1 + varTerm2
  out <- c(CPE, CPEsmooth, sqrt(varCPE))
  names(out) <- c("CPE", "smooth.CPE", "se.CPE")
  out
}
