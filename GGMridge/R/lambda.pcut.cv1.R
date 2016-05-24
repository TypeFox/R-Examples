#----------------------------------------------------------------------#
# Chooses the tuning parameter of the ridge inverse and thresholding   #
# level of the empirical p-values.                                     #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#  train  An n x p data matrix from which the model is fitted.         #
#                                                                      #
#  test   An m x p data matrix from which the model is evaluated.      #
#                                                                      #
#  lambda A vector of candidate tuning parameters.                     #
#                                                                      #
#  pcut   A vector of candidate cutoffs of pvalues.                    #
#                                                                      #
#  Outputs :                                                           #
#                                                                      #
# The total prediction error for all of the candidate lambda           #
#  pvalue cutoff values.                                               #
#----------------------------------------------------------------------#
lambda.pcut.cv1 <- function(train, test, lambda, pcut) {

  p <- ncol(test)
  dp <- diag(p)

  k <- length(lambda)
  lpc <- length(pcut)

  w.upper <- which(upper.tri(dp))
  w.array <- which(upper.tri(dp), arr.ind = TRUE)
  lw <- length(w.upper)

  tunings <- matrix(data = NA, ncol = 2L, nrow = k * lpc)

  S <- cor(train)

  R <- matrix(data = NA, nrow = lw, ncol = k)

  kk <- 0L
  for( la in lambda ) {

    kk <- kk + 1L

    tmp <- try(solve(S + la * dp), silent = TRUE)
    if( is(tmp, "try-error") ) {
      stop("Unable to invert matrix.", call. = FALSE)
    }

    R[,kk] <- -scaledMat(x = tmp)[w.upper]

  }

  transR <- transFisher(x = R)

  efron <- sapply(X = 1L:k,
                  FUN = function(i){getEfronp(z = transR[,i])$correctp})

  risk <- matrix(data = NA, 
                 nrow = k, 
                 ncol = lpc, 
                 dimnames = list(lambda, pcut))

  for( i in 1L:lpc ) {

    thR <- R * {efron < pcut[i]}

    coef.mat <- matrix(data = NA, nrow = p*p, ncol = k)

    for(kk in 1L:k) {

      notzero <- abs(thR[,kk]) > 1.5e-8

      E <- w.array[notzero]
      fit <- structuredEstimate(x = train, E = E)     
      temp <- sqrt(diag(fit$K))
      coef.mat[,kk] <- diag(temp) %*% fit$R %*% diag(1.0/temp)
    }
   
    CV <- 0.0
    for( j in 1L:p ) {
      jp <- j * p
      y <- test[,j]
      D <- test[,-j,drop=FALSE]
      loss <- sapply(X = 1L:k, 
                     FUN = function(i, y, D, coef) {
                             (y - D %*% (coef[,i])[-1L])^2
                           },
                     y = y,
                     D = D,
                     coef = coef.mat[{jp - p + 1L}:jp,,drop=FALSE])

      CV <- CV + colSums(loss) 
    }

    risk[,i] <- CV / length(test) 
  }

  return(risk)
}

