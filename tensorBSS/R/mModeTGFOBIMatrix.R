mModeTGFOBIMatrix <-
function(x, m, lag){
  xm <- mFlatten(x, m)
  mTGFOBIMatrix(xm, lag)
}
