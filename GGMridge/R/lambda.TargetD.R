#----------------------------------------------------------------------#
# Estimation of a weighted average of a sample covariance (correlation)#
# matrix and an identity matrix.                                       #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#  x  Centered data for covariance shrinkage and standardized data for #
#     correlation shrinkage.                                           #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#  The estimates of shrinkage intensity.                               #
#----------------------------------------------------------------------#
lambda.TargetD <- function(x) {

  n <- nrow(x)

  wmean2 <- {t(x) %*% x}^2

  x2 <- x^2
  x2 <- t(x2) %*% x2

  diag.w <- which( row(wmean2) == col(wmean2) )

  varw <- sum(({ x2  - wmean2 / n} * n / {{n-1}^3})[-diag.w])

  S2 <- sum( wmean2[-diag.w] / {n - 1L}^2 )

  lambda <- max(0.0, min( varw / S2, 1.0))

  return(lambda)

}
