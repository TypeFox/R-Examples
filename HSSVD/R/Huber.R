#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Huber : obtain matrix of Huber functions used in initialization              #
#   algorithm of FIT-SSVD, Algorithm 2 arXiv:1112.2433                         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x    : matrix of observed data                                             #
#   beta : degree of "Huberization"                                            #
# Outputs                                                                      #
#   Matrix of Huber functions                                                  #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
Huber <- function(x,
                  beta){

  X.abs <- abs(x)

  #--------------------------------------------------------------------------#
  # delta : the beta-quantile of the absolute values of all the entries in X #
  #--------------------------------------------------------------------------#
  delta <- quantile(x=X.abs, probs=beta)

  #--------------------------------------------------------------------------#
  # Create Huber rho function                                                #
  #--------------------------------------------------------------------------#
  Y <- matrix(0,nrow=nrow(x),ncol=ncol(x))

  matrix.le <- X.abs <= delta

  Y[matrix.le] <- (X.abs[matrix.le])^2
  Y[!matrix.le] <- (2*X.abs[!matrix.le]*delta - delta^2)

  return(Y)
}
