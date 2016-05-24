cvcompute = function(mat, foldid, nlams) {
  ## computes the weighted mean and SD within folds, and
  ##   hence the standard error of the mean
  nfolds = max(foldid)
  outmat = matrix(NA, nfolds, ncol(mat))
  good = matrix(0, nfolds, ncol(mat))
  mat[is.infinite(mat)] = NA
  for (i in seq(nfolds)) {
    mati = mat[foldid == i, ]
    outmat[i, ] = colMeans(mati, na.rm=TRUE)
    good[i, seq(nlams[i])] = 1
  }
  N = colSums(good)
  list(cvraw = outmat, N = N)
}

dwdloss = function(u, qval) {
  ## generalized DWD Loss
  ifelse(u > (qval/(qval+1)), 
    1/(u^qval) * (qval^qval) / ((qval+1)^(qval+1)), 1 - u )
}

err = function(n, maxit) {
  if (n == 0) msg = ""
  if (n < 0) {    
    msg = paste0("convergence for ", -n, 
      "th lambda value not reached after maxit=", maxit, 
      " iterations; solutions for larger lambdas returned")
    n = -1
    msg = paste("From kerneltool fortran code:", msg)
  }
  list(n = n, msg = msg)
}

error.bars = function(x, upper, lower, width=0.02, ...) {
  xlim = range(x)
  barw = diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

getmin = function(lambda, cvm, cvsd) {
  cvmin = min(cvm, na.rm=TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm=TRUE)
  cvmin2 = min(cvm[!is.na(cvsd)])
  lambda.min2 = max(lambda[cvm[!is.na(cvsd)] <= cvmin2], na.rm=TRUE)
  idmin = match(lambda.min2, lambda)
  semin = (cvm + cvsd)[idmin]
  idmin = cvm[!is.na(cvsd)] <= semin
  lambda.1se = max(lambda[idmin])
  id1se = match(lambda.1se, lambda)
  cv.1se = cvm[id1se]
  list(lambda.min = lambda.min, lambda.1se = lambda.1se, 
	  cvm.min = cvmin, cvm.1se = cv.1se)
}

lumloss = function(u, aval, cval) {
  ## LUM Loss
  ifelse(u < (cval/(cval+1)), 1 - u,
    (aval / ((1 + cval) * u - cval + aval) ^ aval ) / (1 + cval))
}
