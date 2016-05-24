#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# rank_est  Perform Wold- or Gabriel-style cross-validation for determining    #
#  the appropriate rank for SVD approximation of a matrix.                     #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#  x      : matrix whose rank is to be estimated                               #
#  method : method of cross-validation (Default: 'wold')                       #
# Outputs                                                                      #
#  estimated rank of x                                                         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
rank_est <- function(x, 
                     method="wold"){

  if( method == "gb" ){
    rankest <- cv.svd.gabriel(x=x, 
                 krow = min(10,floor(nrow(x)/2)), 
                 kcol = min(10,floor(ncol(x)/2)),
                 maxrank = floor(min(nrow(x), ncol(x))/2))
  } else  if( method == "wold"){
    rankest <-  cv.svd.wold(x=x,  
                 k = 5,  
                 maxrank = min(20,ncol(x),nrow(x)),  
                 tol = 1e-4,  
                 maxiter = 20)
  }

  k1 <- which.min(colMeans(rankest$msep))-1

  if(k1 == 0) k1 <- 1

  return(k1)
}
