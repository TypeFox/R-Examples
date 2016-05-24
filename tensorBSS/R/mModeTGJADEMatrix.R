mModeTGJADEMatrix <-
function(x, m, i, j, lag, pm){
  xm <- mFlatten(x, m)
  matCov <- matrixCovariance(xm)
  Eij <- matrix(0, pm, pm)
  Eij[i, j] <- 1
  mat1 <- mTGJADEMatrix(xm, i, j, c(0, lag, lag, 0))
  mat2 <- mTGJADEMatrix(xm, i, j, c(0, lag, 0, lag))
  mat3 <- mTGJADEMatrix(xm, i, j, c(lag, lag, 0, 0))
  return(mat1 + mat2 - mat3 - matCov%*%(Eij + t(Eij) + dim(xm)[2]*(1*(i==j))*diag(pm))%*%matCov)
}
