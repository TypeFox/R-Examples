#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# initialization : Algorithm 2 of FIT-SSVD arXiv:1112.2433                     #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x      : Observed data matrix                                              #
#   beta   : Degree of "Huberization" [0,1]                                    #
#   alpha  : Significance level of a selection test                            #
#   method : Method of cross-validation (Default="wold"                        #
# Outputs                                                                      #
#    u0       : matrix containing the left singular vectors                    #
#    v0       : matrix containing the right singular vectors                   #
#    r_est    : the estimated rank                                             #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
initialization <- function(x, 
                           beta=0.95, 
                           alpha=0.05, 
                           method="wold"){

  #--------------------------------------------------------------------------#
  # Calculate matrix of Huber rho functions                                  #
  #--------------------------------------------------------------------------#
  Y <- Huber(x=x, beta=beta)

  #--------------------------------------------------------------------------#
  # Identify "significant rows"                                              #
  #--------------------------------------------------------------------------#
  index_I <- select.indices(t.vector=rowSums(Y), alpha=alpha)

  #--------------------------------------------------------------------------#
  # Identify "significant columns"                                           #
  #--------------------------------------------------------------------------#
  index_J <- select.indices(t.vector=colSums(Y), alpha=alpha)

  #--------------------------------------------------------------------------#
  # Submatrix includes only significant components                           # 
  #--------------------------------------------------------------------------#
  X_dense <- matrix(x[index_I,index_J],ncol=length(index_J))

  #--------------------------------------------------------------------------#
  # Obtain rank and SVD of X_dense                                           #
  #--------------------------------------------------------------------------#
  if( (ncol(X_dense) > 1) && (nrow(X_dense) > 1) ){
    r_est <- rank_est(x=X_dense, method=method)
  } else {
    r_est <- 1
    X_dense <- matrix(X_dense, ncol=1)
  }
  U0 <-  matrix(0,nrow=nrow(x),ncol=r_est)
  V0 <-  matrix(0,nrow=ncol(x),ncol=r_est)

  svd_dense <- svd(x=X_dense)

  U0[index_I,] <- svd_dense$u[,1:r_est]
  V0[index_J,] <- svd_dense$v[,1:r_est]

  return(list(      u0 = U0,
                    v0 = V0,
                 r_est = r_est))
}
