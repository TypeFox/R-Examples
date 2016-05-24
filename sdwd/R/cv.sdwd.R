cv.sdwd = function(x, y, lambda=NULL, nfolds=5, foldid, ...) {
  N = nrow(x)
  ## Fit the model once to get dimensions etc of output
  y = drop(y)
  sdwd.object = sdwd(x, y, lambda=lambda, ...)
  lambda = sdwd.object$lambda
  # predict -> coef
  nz = sapply(coef(sdwd.object, type="nonzero"), length)
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), 
      length=N)) else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  ## Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which = foldid == i
    outlist[[i]] = sdwd(x=x[!which, , drop=FALSE], 
                        y=y[!which], lambda=lambda, ...)
  }
  ## What to do depends on the model fit
  pred.loss = "misclass"
  fun = paste("cv", class(sdwd.object)[[2]], sep=".")
  cvstuff = do.call(fun, list(outlist, lambda, x, y, foldid))
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  cvname = cvstuff$name
  out = list(lambda=lambda, cvm=cvm, cvsd=cvsd, 
              cvupper=cvm+cvsd, cvlower=cvm - cvsd, nzero=nz, 
              name=cvname, sdwd.fit=sdwd.object)
  lamin = getmin(lambda, cvm, cvsd)
  obj = c(out, as.list(lamin))
  class(obj) = "cv.sdwd"
  obj
} 

cv.sdwdNET = function(outlist, lambda, x, y, foldid) {
  ### Turn y into c(-1,1)
  y = as.factor(y)
  y = c(-1, 1)[as.numeric(y)]
  nfolds = max(foldid)
  predmat = matrix(NA, length(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    preds = predict(fitobj, x[which, , drop=FALSE], type="link")
    nlami = length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  pred.loss="misclass"
  typenames = c(misclass = "Mis-classification Error")
  N = length(y) - apply(is.na(predmat),2,sum)
  cvraw = switch(pred.loss, 
    misclass = (y != ifelse(predmat > 0, 1, -1)))
  if (length(y)/nfolds >= 3) {
    cvob = cvcompute(cvraw, foldid, nlams)
    cvraw = cvob$cvraw; N = cvob$N
  }     
  cvm = apply(cvraw, 2, mean, na.rm=TRUE)
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, 
    mean, na.rm=TRUE)/(N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 

predict.cv.sdwd = function(object, newx, s=c("lambda.1se", 
    "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda = s else if (is.character(s)) {
      s = match.arg(s)
      lambda = object[[s]]
    } else stop("Invalid form for s")
  predict(object$sdwd.fit, newx, s=lambda, ...)
} 


plot.cv.sdwd = function(x, sign.lambda=1, ...) {
  cvobj = x
  xlab = "log(Lambda)"
  if (sign.lambda < 0) 
    xlab = paste0("-", xlab)
  plot.args = list(x = sign.lambda * log(cvobj$lambda), 
    y = cvobj$cvm, type = "n", xlab = xlab, ylab = cvobj$name,
    ylim = range(cvobj$cvupper, cvobj$cvlo))
  new.args = list(...)
  if (length(new.args)) 
    plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper, 
    cvobj$cvlo, width=0.01, col="darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, 
    pch=20, col="red")
  axis(side=3, at=sign.lambda * log(cvobj$lambda), 
    labels=paste(cvobj$nz), tick=FALSE, line=0)
  abline(v=sign.lambda * log(cvobj$lambda.min), lty=3)
  abline(v=sign.lambda * log(cvobj$lambda.1se), lty=3)
  invisible()
} 

coef.cv.sdwd = function(object, s=c("lambda.1se", 
    "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda = s else if (is.character(s)) {
      s = match.arg(s)
      lambda = object[[s]]
    } else stop("Invalid form for s.")
  coef(object$sdwd.fit, s=lambda, ...)
} 
