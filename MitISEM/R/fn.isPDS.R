# This function tests if the input is a positive definite symmetric matrix
#
# inputs:
#     Sigma   : [kxk matrix] to test pds
# outputs: 
#     r       : [logical], 'TRUE' if Sigma is pds, 'FALSE' otherwise
#
# note: Function tests the stability checking 
#       determinant, Choleski decomposition and eigenvalues
# author : Nalan Basturk
# date   : 20120912

fn.isPDS <- function(Sigma){
  r   <- !all(eigen(Sigma)$values>0)         # not pds matrix?
  if(!r){
    r <- try(chol(Sigma),silent=TRUE)        # issue in Choleski decomp.?
    r = !is.matrix(r)  
  }  
  tmp <- abs(det(Sigma))
  tol <- 1e25
  if(!r)
    r  <-  (tmp>=tol | tmp<=1/tol)           # determinant too small/large?
  if(!r){                                    # issue in matrix inversion?
    r  <- try(solve(Sigma),silent=TRUE)
    r   = !is.matrix(r)
  }
  as.logical(r)
}  
