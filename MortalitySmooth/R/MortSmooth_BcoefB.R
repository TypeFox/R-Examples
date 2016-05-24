MortSmooth_BcoefB <-
function(X1, X2, mat){
  ## Input:
  ## X1: first marginal matrix
  ## X2: second marginal matrix
  ## mat: matrix to be multiplied
  ## Output:
  ## BcoefB: a matrix that multiplies X1 by mat
  ##         and the transpose of X2
    BcoefB <- X1%*%mat%*%t(X2)
    return(BcoefB)
}
