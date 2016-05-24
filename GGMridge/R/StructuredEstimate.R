#----------------------------------------------------------------------#
# Estimation of partial correlation matrix given zero structure.       #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#  x  n x p data matrix; p = number of variables p, n = sample size    #
#                                                                      #
#  E  The row and column indices of zero entries of the partial        #
#     correlation matrix.                                              #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#  R   The partial correlation matrix.                                 #
#                                                                      #
#  K   The inverse covariance matrix.                                  #
#                                                                      #
#  RSS The residual sum of squares.                                    #
#                                                                      #
#----------------------------------------------------------------------#
structuredEstimate <- function(x, E) {

  E <- matrix(E, ncol = 2L)

  n <- nrow(x)
  p <- ncol(x)

  RSS <- numeric(p)

  K <- matrix(data = 0.0, nrow = p, ncol = p)
  diag(K) <- apply(X = x, MARGIN = 2L, FUN = var)

  for( i in 1L:p ) {

    y <- x[,i] 

    ne <- c( E[{E[,1L] == i},2L], E[{E[,2L] == i},1L])
    lne <- length(ne)

    if( lne == 0L ) {
      RSS[i] <- y %*% y
    } else {
      D <- x[, ne, drop = FALSE] 
      svdD <- svd(x = D)
      coef <- svdD$v %*% { t(svdD$v) / svdD$d^2 } %*% t(D) %*% y
      RSS[i] <- sum( {y - D %*% coef}^2)
      resid.var <- {n - lne} / RSS[i]
      coef.cl <- c( resid.var, - coef * resid.var)
      id.cl <- c(i, ne)  
      o.cl <- order(id.cl)
      K[i,id.cl[o.cl]] <- coef.cl[o.cl]
    }
  }

  absK <- abs(K)
  tmp1 <- sqrt(absK * t(absK))
  tmp2 <- sign(K) * upper.tri(K)
  tmp3 <- {tmp2 + t(tmp2)} * tmp1
  diag(tmp3) <- diag(K)
  R <- - scaledMat(x = tmp3)
  diag(R) <- 1.0

  return(list("R" = R, 
              "K" = tmp3,  
              "RSS" = RSS))
}
