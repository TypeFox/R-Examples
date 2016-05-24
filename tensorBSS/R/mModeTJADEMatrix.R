mModeTJADEMatrix <-
function(x, m, i, j){
  xm <- mFlatten(x, m)
  matCov <- matrixCovariance(xm)
  mJADEMatrix(xm, i, j, matCov)
}
