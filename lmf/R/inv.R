inv <-
function (a)
{
  if(!any(is.na(a)))
    #Find the upper triangular matrix of X and estimate the inverse
    #of X from this matrix
    chol2inv(chol(a))
  else
    #If NA values in matrix, return NA matrix
    matrix(NA, nrow = dim(a)[1], ncol = dim(a)[2])
}
