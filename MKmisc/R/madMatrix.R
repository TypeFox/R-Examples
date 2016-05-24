################################################################################
## computes MAD for a matrix - similar to covariance matrix
################################################################################
madMatrix <- function(x){
  anz <- ncol(x)
  mad.matrix <- diag(rep(0, anz))
  for(i in 1:(anz-1)){
    mad.matrix[i, ((i+1):anz)] <- apply(as.matrix(x[,i] - x[,(i+1):anz]), 2, mad, na.rm = TRUE)
    mad.matrix[(i+1):anz, i] <- mad.matrix[i, ((i+1):anz)]
  }
  colnames(mad.matrix) <- colnames(x)

  return(mad.matrix)
}
