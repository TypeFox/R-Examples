#----------------------------------------------------------------------#
# The partial correlation matrix is estimated by p separate ridge      #
# regressions with the parameters selected by cross validation.        #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#  x       n x p data matrix; n = sample size, p = # of variables.     #
#                                                                      #
#  fold    fold-cross validation                                       #
#                                                                      #
#  lambda  the candidate ridge parameters                              #
#                                                                      #
#  verbose TRUE indicates extra printing                               #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#  R          The partial correlation matrix.                          #
#                                                                      #
#  lambda.sel The selected tuning parameters for p ridge regressions.  #
#                                                                      #
#----------------------------------------------------------------------#
R.separate.ridge <- function(x, fold, lambda, verbose = FALSE) {

  if( !is(x, "matrix") ) {
    stop("x must be an object of class matrix.", call. = FALSE)
  }
  if( nrow(x) <= 1L ) {
    stop("More than one sample is required in x.", call. = FALSE)
  }
  if( ncol(x) <= 1L ) {
    stop("More than one covariate is required in x.", call. = FALSE)
  }

  n <- nrow(x)
  p <- ncol(x)
  x <- scale(x = x, center = TRUE, scale = TRUE)
    
  coefNI <- matrix(data = 0.0, nrow = p, ncol = p)
  diag.coefNI <- numeric(p)
  lambda.sel <- numeric(p)

  for( i in 1L:p ) {

    if( verbose ) cat("variable=", i, " ", date(), "\n")

    tempX <- x[,-i,drop = FALSE]
    tempY <- x[,i]

    fit.lambda <- ne.lambda.cv(y = tempY,
                               x = tempX,
                               lambda = lambda,
                               fold = fold) 

    lambdai <- fit.lambda$lambda[which.min(fit.lambda$spe)]
    lambda.sel[i] <- lambdai

    ridgefit <- try(lm.ridge(tempY ~ tempX - 1, lambda = lambdai), 
                    silent = TRUE)
    if( is(ridgefit, "try-error") ) {
      stop("Ridge regressioni not successful.", call. = FALSE)
    }

    predY <- tempX %*% coef(ridgefit) 
    Ds <- try(svd(x = tempX), silent = TRUE)
    if( is(Ds, "try-error") ) {
      stop("Unable to obtain singular value decomposition.", 
           call. = FALSE)
    }
    dvals2 <- Ds$d * Ds$d

    diag.coefNI[i] <- {n - sum(dvals2 / {dvals2 + lambdai})} / 
                      sum({tempY - predY}^2)
    coefNI[i,-i] <- -diag.coefNI[i] * coef(ridgefit)
  }

  absc <- abs(coefNI)
  tmp1 <- sqrt(absc * t(absc))
  tmp2 <- sign(coefNI) * upper.tri(coefNI)
  tmp3 <- (tmp2 + t(tmp2)) * tmp1
  diag(tmp3) <- diag.coefNI
  R <- - scaledMat(x = tmp3)
  diag(R) <- 1.0

  return(list("R" = R,
              "lambda.sel" = lambda.sel))
}
