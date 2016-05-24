#####################################################################
#  This is the general calling program for fitting GLARMA models with
#  various distributions (Poisson, Binomial, Negative Binomial),
#  method of scoring or Newton Raphson, Pearson Residuals or Score
#  residuals.
#
#
#   User Manual:
#   Dunsmuir (2010) "R software for fitting observation driven
#   regression models for univariate time series"
#
#   Contributers: See Dunsmuir(2010).
#
#   Last Modified: March 29, 2011.
#
####################################################################



glarma <- function(y, X, offset = NULL, type = "Poi", method = "FS",
                   residuals = "Pearson",
                   phiLags,  thetaLags, phiInit, thetaInit, beta, alphaInit,
                   alpha = 1, maxit = 30, grad = 2.22e-16) {

  call <- match.call()
  phi <- phiGen(phiLags, phiInit)
  phiLags <- phi[[1]]
  phi.init <- phi[[2]]

  theta <- thetaGen(thetaLags, thetaInit)
  thetaLags <- theta[[1]]
  theta.init <- theta[[2]]

  delta <- deltaGen(y = y, X = X, offset = offset,
                    phiInit = phi.init, thetaInit = theta.init,
                    type = type, alpha = alpha, beta = beta,
                    alphaInit = alphaInit)


  deltanames <- names(delta)
  min        <- 1
  iter       <- 0
  ErrCode    <- 0
  Werror     <- 0

  while ((iter < maxit) & (min > grad) & (ErrCode == 0)) {
    if (type == "Poi" & residuals == "Pearson")

      GL <- glarmaPoissonPearson(y, X, offset = offset, delta,
                                 phiLags, thetaLags,
                                 method = method)

    if (type == "NegBin" & residuals == "Pearson")
      GL <- glarmaNegBinPearson(y, X, offset = offset, delta,
                                phiLags, thetaLags,
                                method = method)

    if (type == "Bin" & residuals == "Pearson")
      GL <- glarmaBinomialPearson(y, X, offset = offset, delta,
                                  phiLags, thetaLags,
                                  method = method)

    if (type == "Poi" & residuals == "Score")
      GL <- glarmaPoissonScore(y, X, offset = offset, delta,
                               phiLags, thetaLags,
                               method = method)

    if (type == "NegBin" & residuals == "Score")
      GL <- glarmaNegBinScore(y, X, offset = offset, delta,
                              phiLags, thetaLags,
                              method = method)

    if (type == "Bin" & residuals == "Score")
      GL <- glarmaBinomialScore(y, X, offset = offset, delta,
                                phiLags, thetaLags,
                                method = method)

    if (type == "Bin" & residuals == "Identity")
      GL <- glarmaBinomialIdentity(y, X, offset = offset, delta,
                                   phiLags, thetaLags,
                                   method = method)

    temp <- mySolve(-GL$ll.dd)

    if (temp$ErrCode == 0) {
      GL$cov    <- temp$Ainv
      step      <- GL$cov %*% GL$ll.d
      delta.old <- GL$delta
      delta     <- delta.old + step
      min       <- max(abs(GL$ll.d))
      iter      <- iter + 1
    }
    ErrCode <- temp$ErrCode
    if (ErrCode == 1) {
        stop(paste(switch(method, "FS" = "Fisher Scoring",
                          "NR" = "Newton Raphson"),
                   "fails to converge from the initial estimates."))
    }
  }


  ## keep model structure
  names(GL)[which(names(GL) == "ll")] <- "logLik"
  names(GL)[which(names(GL) == "ll.d")] <- "logLikDeriv"
  names(GL)[which(names(GL) == "ll.dd")] <- "logLikDeriv2"
  names(GL)[which(names(GL) == "e")] <- "residuals"
  GL$null.deviance <- initial(y, X, offset = offset, type = type,
                              alpha = alpha)$null.deviance
  GL$df.null <- initial(y, X, offset = offset, type = type,
                        alpha = alpha)$df.null
  GL$r         <- ncol(X)
  GL$pq        <- length(phiLags) + length(thetaLags)
  GL$phiLags   <- phiLags
  GL$thetaLags <- thetaLags
  GL$y         <- y
  GL$X         <- X
  GL$offset    <- offset
  GL$type      <- type
  GL$method    <- method
  GL$residType <- residuals
  GL$call      <- call
  ## computes the results
  names(delta) <- deltanames
  GL$delta     <- delta
  GL$iter      <- iter
  GL$errCode   <- ErrCode
  if (min(is.finite(GL$W)) == 0) Werror <- 1
  GL$WError <- Werror
  GL$min    <- min
  GL$aic    <- -2 * GL$logLik + 2 * (GL$pq + GL$r)

  class(GL)  <- "glarma"
  GL
}

## print method for glarma function
print.glarma <- function(x, ...) {
   cat("\nCall: ", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")
   if (x$type == "NegBin"){
     cat("Negative Binomial Parameter:\n")
     print.default(format(x$delta[length(x$delta)]),
                   print.gap = 2, quote = FALSE)
     if (x$pq > 0){
       cat("\nGLARMA Coefficients:\n")
       print.default(format(x$delta[(ncol(x$X) + 1):(ncol(x$X) +
                                                     length(x$thetaLags) +
                                                     length(x$phiLags))]),
                     print.gap = 2, quote = FALSE)
     }
   } else{
     if (x$pq > 0){
         cat("GLARMA Coefficients:\n")
         print.default(format(x$delta[(ncol(x$X) + 1):length(x$delta)]),
                       print.gap = 2, quote = FALSE)
     }
   }
   cat("\nLinear Model Coefficients:\n")
   print.default(format(x$delta[1:ncol(x$X)]), print.gap = 2,
                 quote = FALSE)
   cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
       NROW(x$y) - NROW(x$delta), "Residual")
   cat("\nNull Deviance:", x$null.deviance,
       "\nResidual Deviance:", sum(x$residuals^2),
       "\nAIC:", format(x$aic), "\n\n")
}

