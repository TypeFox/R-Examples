cv.kerndwd = function(x, y, kern, lambda, nfolds=5, foldid, ...) {
  ####################################################################
  ## data setup
  y = drop(y)
  x = as.matrix(x)
  x.row = as.integer(NROW(x))
  if (length(y) != x.row) 
    stop("x and y have different number of observations.")  
  ####################################################################
  ## fit DWD nfold times and save them 
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), length=x.row)) else 
    nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be larger than 3; nfolds=5 recommended.")
  lambda = sort(lambda, decreasing=TRUE)
  nlams = double(nfolds)
  outlist = as.list(nlams)
  for (i in seq(nfolds)) {
    which = foldid == i
    outlist[[i]] = kerndwd(x=x[!which, , drop=FALSE], y=y[!which], 
      kern=kern, lambda=lambda, ...)
    if (outlist[[i]]$jerr != 0) 
      stop(paste("Error occurs when fitting the", i, "th folder."))
  }
  cvstuff = cvpath(outlist, x, y, kern, lambda, foldid, 
    x.row, ...)
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  ## wrap up output
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
    cvup = cvm + cvsd, cvlo = cvm - cvsd, 
    name = cvstuff$name)
  obj = c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) = "cv.kerndwd"
  obj
} 

cvpath = function(outlist, x, y, kern, lambda,  
  foldid, x.row, ...){
  ####################################################################
  if (length(levels(factor(y))) == 2)
    y = c(-1, 1)[as.factor(y)]
  # if (!all(y %in% c(-1, 1))) 
  #   stop("y should be a factor with two levels.")
  typenames = "mis-classification error"
  predmat = matrix(NA, x.row, length(lambda))
  nfolds = max(foldid)
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    preds = predict(fitobj, kern, x[!which, , drop=FALSE], 
      x[which, , drop=FALSE], type="link")
    nlami = length(fitobj$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  ## select the lambda according to predmat
  cvraw = (y != ifelse(predmat > 0, 1, -1))
  if (length(y)/nfolds >= 3) {
    cvob = cvcompute(cvraw, foldid, nlams)
    cvraw = cvob$cvraw
    cvn = cvob$N
  } else cvn = length(y) - colSums(is.na(predmat))    
  cvm = colMeans(cvraw, na.rm=TRUE)
  cvsd = sqrt(colMeans(scale(cvraw, cvm, FALSE)^2, 
    na.rm=TRUE)/(cvn - 1))
  out = list(cvm = cvm, cvsd = cvsd, 
    name = "mis-classification error")
  out
}
