print.artfima <-
function(x, ...){
  glp <- x$glp
  p <- x$arimaOrder[1]
  d <- x$arimaOrder[2]
  q <- x$arimaOrder[3]
  est <- numeric(0)
  estT <- character(0)
  glpOrder <- x$glpOrder
  est <- c(est, x$constant)
  estT <- c(estT, paste0(
            ifelse(d==0,"mean","constant"),
            ifelse(x$blueQ," (blue)", "")
            ))
  if (glpOrder==2) {
    est <- c(est, x$lambdaHat)
    estT <- c(estT, "lambda")
    est <- c(est, x$dHat)
    estT <- c(estT, "d")
  } else {
    if(glpOrder==1) {
      est <- c(est, x$dHat)
      estT <- c(estT, "d")
    }
  }
  if(p>0) {
    est <- c(est, x$phiHat)
    estT <-c(estT, paste0("phi(", paste0(1:p, ")")))
  }
  if(q>0) {
    est <- c(est, x$thetaHat)
    estT <- c(estT, paste0("theta(", paste0(1:q, ")")))      
  }
  whichModel <- paste0(x$glp, "(", p, ",", d, ",", q, ")" )
  cat(paste0(whichModel, ", MLE Algorithm: ", x$likAlg),fill=TRUE)
  cat(paste0("snr = ", round(x$snr,3), ", sigmaSq = ", 
             x$sigmaSq), fill=TRUE)
  if (x$convergence!=0) {
    cat(paste0("Note: possible problem with convergence\n  convergence=",
               x$convergence), fill=TRUE)
    cat(paste("Algorithm used: ", x$algorithm), fill=TRUE)
    if (!is.null(x$message)) 
      cat(paste0("message =", x$message), fill=TRUE)
  }
  k <- length(est)
  LL <- x$LL
  aic <- -2*LL + 2*k
  bic <- -2*LL + k*log(x$n)
  cat(paste0("log-likelihood =", round(LL,3), ", AIC = ", round(aic,1), 
      ", BIC = ", round(bic,1)),  fill=TRUE)
  if (x$onBoundary) {
    cat("Warning: estimates converged to boundary!", fill=TRUE)
  }
  print(matrix(c(est, x$seMean, x$se), ncol=2, 
         dimnames=list(estT, c("est.", "se(est.)"))))   
}
