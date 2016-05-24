orderedLasso = function(x, y, lambda, intercept = TRUE, b0 = NULL, beta_pos = NULL, beta_neg = NULL, 
                        method = c("Solve.QP", "GG"), strongly.ordered = FALSE,
                        standardize = TRUE, trace = FALSE, niter = 500, iter.gg = 100,  epsilon = 1e-8){  
  
  stopifnot(nrow(x) == length(y), lambda >= 0)
  stopifnot(class(lambda) == "numeric")
  stopifnot(is.finite(x), is.finite(y), is.finite(lambda))
  if (is.null(beta_pos)) beta_pos = rep(0, ncol(x))
  if (is.null(beta_neg)) beta_neg = rep(0, ncol(x))
  if(missing(method)) method = "Solve.QP"
  this.call <- match.call()

  if (standardize) {
    stdeviation_x =  apply(x, 2, sd) 
    stdeviation_inverse = 1 / stdeviation_x   
  }
  if (intercept){
    mean_y = mean(y)
    y = y - mean_y
    mean_x = apply(x, 2, mean)
  }
  x = scale(x, intercept, standardize)
  xx = -x
  if(trace){
    val =   ObjMinFun(x = x, y = y, bp = beta_pos, bn = beta_neg, lambda = lambda)
    cat(c("initial", val), fill = T)  
  }

  if(method == "GG"){
    ii = 0
    go = TRUE
    while(go & ii< niter){
      go=FALSE
      ii = ii + 1 
      beta_pos_old = beta_pos
      beta_neg_old = beta_neg
      r = y - xx %*% beta_neg
      junk = fastgg(x, r, +1, beta_pos, beta_neg, lam  = lambda,  trace = trace, inneriter.gg = iter.gg)
      beta_pos= junk$beta_pos
      beta_neg= junk$beta_neg 
      if(trace){
        val =  ObjMinFun(x = x, y = y, bp = beta_pos, bn = beta_neg, lambda = lambda)
        cat(c("aft pos",val), fill=T)
      }
      r = y - x %*% beta_pos  
      junk = fastgg(xx, r, -1,beta_pos, beta_neg, lam = lambda, trace=trace, inneriter.gg = iter.gg)  
      
      beta_pos = junk$beta_pos
      beta_neg = junk$beta_neg
      
      if(trace){
        val =  ObjMinFun(x = x, y = y, bp = beta_pos, bn = beta_neg, lambda = lambda)
        cat(c("aft neg",val), fill=T)
      }
      if (sum(abs(beta_pos_old-beta_pos))>epsilon | sum(abs(beta_neg_old - beta_neg)) > epsilon) go=TRUE
    }
  } 
  if(method =="Solve.QP"){
    update = ordLas2(x, y, lam = lambda)
    beta_pos = update$bp
    beta_neg = update$bn
    val = ObjMinFun(x = x, y = y, bp = beta_pos, bn = beta_neg, lambda = lambda)
    if(trace) cat("after", val, fill = T)
  }
  beta = beta_pos - beta_neg
  ordered.test = FALSE
  beta.diff = beta[1:(length(beta)-1)]-beta[2:(length(beta))]
  if(any(beta.diff < (-1e-5))) ordered.test = TRUE
  if(strongly.ordered & ordered.test){
    signvec = sign(beta)
    beta.ordered = ordLasSignPos(x = x,y =y, lam = lambda, signvec = signvec)$bp
  }
  if(strongly.ordered == TRUE & (!ordered.test)) beta.ordered = beta
  if(!strongly.ordered) beta.ordered = b0.ordered =fitted.ordered =  NULL
  if (intercept){
    if (standardize){
      beta_pos = stdeviation_inverse * beta_pos
      beta_neg = stdeviation_inverse * beta_neg
      beta = beta_pos - beta_neg
      b0 = as.numeric(mean_y - mean_x %*% beta)   
      if(strongly.ordered){
        beta.ordered = stdeviation_inverse * beta.ordered
        b0.ordered = as.numeric(mean_y - mean_x %*% beta.ordered)
      }
    }else{
      b0 = as.numeric(mean_y - mean_x %*% beta)
      if(strongly.ordered){
        b0.ordered = as.numeric(mean_y - mean_x %*% beta.ordered)
      }
    }
    fitted = x %*% beta + mean_y
    if(strongly.ordered) fitted.ordered = x %*% beta.ordered + mean_y
  }
  if (!intercept){
    if(standardize){
      beta_pos = stdeviation_inverse * beta_pos
      beta_neg = stdeviation_inverse * beta_neg 
      beta = beta_pos - beta_neg
      b0 = NULL
      if(strongly.ordered){
        b0.ordered = NULL
        beta.ordered = stdeviation_inverse * beta.ordered
      }
    }else{
      b0 = NULL
      if(strongly.ordered) b0.ordered = NULL
    }
    fitted = x %*% beta
    if(strongly.ordered) fitted.ordered = x%*% beta.ordered 
  }  
  if(!strongly.ordered) fitted.ordereed = b0.ordered = beta.ordered = NULL
  type = "gaussian"
  out = list(bp = beta_pos, bn = beta_neg, beta = beta, b0 = b0, b0.ordered = b0.ordered, beta.ordered = beta.ordered, 
             strongly.ordered = strongly.ordered, type = type, call = this.call, method = method, fitted = fitted, fitted.ordered = fitted.ordered)
  val =  ObjMinFun(x = x, y = y, bp = beta_pos, bn = beta_neg, lambda = lambda)
  class(out) <- "orderedLasso"
  out
} 
orderedLasso.path = function(x, y, lamlist = NULL, minlam = NULL, maxlam = NULL, nlam = 50,  flmin = 5e-3, intercept = TRUE, 
                             standardize = TRUE, method = c("Solve.QP", "GG"),niter = 500, iter.gg = 100, strongly.ordered  = FALSE, 
                             trace = FALSE, epsilon = 1e-5){     
  stopifnot(nrow(x) == length(y))
  stopifnot(is.finite(x), is.finite(y))
  this.call = match.call()

  if(missing(method)) method = "Solve.QP"
  
  if (is.null(maxlam)){
    if (!is.null(minlam)) stop("Cannot have maxlam = NULL if minlam is non-null.")
    y_temp = y - mean(y)
    maxlam <- max(abs(crossprod(x,y_temp)))
    minlam <- maxlam * flmin
  }
 
  if (is.null(minlam)) minlam <- maxlam * flmin
  if (is.null(lamlist))  lamlist <- exp(seq(log(maxlam),log(minlam),length=nlam))
 
  nlam <- length(lamlist)
  p = ncol(x)
  beta_pos_lambda = matrix(NA, ncol = nlam, nrow = p)
  beta_neg_lambda = matrix(NA, ncol = nlam, nrow = p)
  beta = matrix(NA, ncol =  nlam, nrow = p) 

  if(strongly.ordered) beta.ordered = matrix(NA, ncol = nlam, nrow = p)
  fitted  = matrix(NA, ncol =  nlam, nrow = nrow(x))
  if(strongly.ordered) fitted.ordered = matrix(NA, ncol =  nlam, nrow = nrow(x))
  bpinit = bninit=rep(NA, p) 
  err = rep(NA, length =  nlam)

  if(strongly.ordered) err.ordered = rep(NA, length =  nlam)

  if(intercept) {
    b0 = rep(NA, length =  nlam) 
    b0init = 0
    if(strongly.ordered) b0.ordered = rep(NA, length = nlam)
  }else{
    b0 = b0init = b0.ordered = NULL
  }
  
  for (i in (1: nlam)){
    if (trace){
       cat(c("lambda=",lamlist[i]),fill=T)
    }
    if(i ==1){
      bpinit = rep(0, length = p)
      bninit = rep(0, length = p)
      if(intercept) b0init = 0
    }

    if(i>1) {
      bpinit = beta_pos_lambda[, i-1]
      bninit = beta_neg_lambda[, i-1]
      if(intercept) b0init = b0[i-1]
    } 
    est = orderedLasso(x = x, y = y, lambda = lamlist[i], intercept = intercept, b0 = b0init, beta_pos = bpinit, beta_neg = bninit, 
                             standardize = standardize, trace= FALSE, niter = niter, iter.gg = iter.gg, epsilon = epsilon, 
                             method = method, strongly.ordered = strongly.ordered)
    beta_pos_lambda[, i] = est$bp
    beta_neg_lambda[, i] = est$bn
    if(strongly.ordered)  beta.ordered[, i] = est$beta.ordered
    
    if(intercept){
      b0[i] = est$b0
      if(strongly.ordered) b0.ordered[i] = est$b0.ordered
    }

    beta[,i] = est$beta
    fitted[, i] = est$fitted

    if(strongly.ordered) fitted.ordered[,i] = est$fitted.ordered
    err[i] = mean((y - fitted[,i])^2)
    if(strongly.ordered) err.ordered[i] = mean((y - fitted.ordered[,i])^2)
  }
  
  if(!strongly.ordered) beta.ordered = fitted.ordered= err.ordered = b0.ordered = NULL

  out = list(bp = beta_pos_lambda, bn = beta_neg_lambda, beta = beta, b0 = b0, b0.ordered = b0.ordered, beta.ordered = beta.ordered, strongly.ordered = strongly.ordered,
            lamlist = lamlist, err = err, err.ordered = err.ordered, method = method, call = this.call)
  class(out) = "orderedLasso.path"
  return(out)
}
 
orderedLasso.cv = function(x, y, lamlist = NULL, minlam = NULL, maxlam = NULL, nlam = 50, 
                          flmin = 5e-4, strongly.ordered = FALSE, 
                          intercept = TRUE, standardize = TRUE, nfolds=10, folds=NULL,  niter = 500, iter.gg = 100, 
                          method = c("Solve.QP", "GG"),trace = FALSE, epsilon = 1e-5){
  length_y <- length(y)
  stopifnot(nrow(x) == length_y)
  stopifnot(is.finite(x), is.finite(y))
  if(missing(method)) method = "Solve.QP"
  this.call = match.call()
  errfun= function(yhat, y){(yhat - y)^2}

 if (is.null(maxlam)) {
    if (!is.null(minlam)) stop("Cannot have maxlam = NULL if minlam is non-null.")
    y_temp = y - mean(y)
    maxlam =  max(abs(crossprod(x,y_temp)))
    minlam <- flmin * maxlam
  }

  if (is.null(minlam)) minlam <- maxlam * flmin
  fit = NULL
  if (is.null(lamlist))  lamlist <- exp(seq(log(minlam), log(maxlam),length=nlam))
  nlam <- length(lamlist)

  if(is.null(folds)) {
    folds <- split(sample(1:length_y), rep(1:nfolds, length = length_y))
  }else {
    stopifnot(class(folds)=="list")
    nfolds <- length(folds)
  } 
  err2 = matrix(NA,nrow=nfolds,ncol=length(lamlist))
  if(strongly.ordered) err2.ordered = matrix(NA,nrow=nfolds,ncol=length(lamlist))
 
  fit = orderedLasso.path(x, y, lamlist = lamlist, minlam = minlam, maxlam = maxlam, 
                          nlam = nlam, flmin = flmin, intercept = intercept, strongly.ordered = strongly.ordered,
                          standardize = standardize, niter = niter, iter.gg = iter.gg, trace = FALSE, epsilon = epsilon, method = method)
 
  nonzeros = colSums((fit$bp-fit$bn) !=0)
  if(strongly.ordered) nonzeros.ordered = colSums((fit$bp-fit$bn)!= 0)

  for (ii in 1: nfolds){
    if(trace){
      cat("Fold", ii, ":\n")
    }
    a = orderedLasso.path(x = x[-folds[[ii]],], y=y[-folds[[ii]]], lamlist = lamlist, 
                                intercept = intercept, standardize = standardize, method = method, strongly.ordered = strongly.ordered,
                                trace = trace, niter = niter, iter.gg = iter.gg, epsilon = epsilon)
    predicted.model = predict.orderedLasso.path(a, newdata = x[folds[[ii]],])

    yhatt = predicted.model$yhat
    if(strongly.ordered) yhatt.ordered = predicted.model$yhat.ordered
    
    temp = matrix(y[folds[[ii]]],nrow=length(folds[[ii]]),ncol=length(lamlist))
    err2[ii, ] = colMeans(errfun(yhatt,temp))
    if(strongly.ordered) err2.ordered[ii, ] = colMeans(errfun(yhatt.ordered,temp))
 
  } 
  errm=colMeans(err2)
  errse=sqrt(apply(err2,2,var)/nfolds)
  o=which.min(errm)
  lamhat=lamlist[o]
  oo = errm >= errm[o] + errse[o]
  lamhat.1se=sort(lamlist[oo & lamlist>=lamhat])[1] 
  if(strongly.ordered){
    errm.ordered=colMeans(err2.ordered)
    errse.ordered=sqrt(apply(err2.ordered,2,var)/nfolds)
    o.ordered=which.min(errm.ordered)
    lamhat.ordered=lamlist[o.ordered]
    oo.ordered = errm.ordered >= errm.ordered[o.ordered] + errse.ordered[o.ordered]
    lamhat.ordered.1se=sort(lamlist[oo.ordered & lamlist>=lamhat.ordered])[1]
  }else{
    errm.ordered = errse.ordered = o.ordered=lamhat.ordered=oo.ordered=lamhat.ordered.1se = nonzeros.ordered= NULL
  }
  obj <- list(lamlist = lamlist, cv.err = errm, cv.err.ordered = errm.ordered, cv.se =errse, cv.se.ordered = errse.ordered, 
              lamhat=lamhat, lamhat.ordered = lamhat.ordered,  strongly.ordered = strongly.ordered,
              lamhat.1se = lamhat.1se, lamhat.ordered.1se = lamhat.ordered.1se, 
              folds=folds, method = method, call = this.call, nonzeros = nonzeros, nonzeros.ordered = nonzeros.ordered)
  class(obj) <- "orderedLasso.cv"
  obj
}

print.orderedLasso.cv = function(x,...){
  cat("Call\n")
  dput(x$call)
  if (x$strongly.ordered){
    mat = cbind(x$lamlist, x$cv.err, x$cv.se, x$cv.err.ordered, x$cv.se.ordered)
    dimnames(mat) = list(NULL, c("Lambda", "Mean CV Error", "SE", "Mean CV Error Ordered", "SE Ordered"))
  }else{
    mat = cbind(x$lamlist, x$cv.err, x$cv.se) 
    dimnames(mat) = list(NULL, c("Lambda", "Mean CV Error", "SE"))
  }
  cat("\n")
  print(mat, quote = FALSE)
  cat("\n")
  if(x$strongly.ordered){ 
    cat(c("lamhat=",round(x$lamhat,5),"lamhat.1se=",round(x$lamhat.1se,5), "lamhat.ordered = ", round(x$lamhat.ordered,5), "lamhat.ordered.1se = ", round(x$lamhat.ordered.1se,5)),fill=TRUE)
  }else{
    cat(c("lamhat=",round(x$lamhat,5),"lamhat.1se=",round(x$lamhat.1se,5)),fill = TRUE)
  }
  mat
}

plot.orderedLasso.path = function(x, ...){
  lambda = x$lamlist
  beta= t(x$beta)
  p = ncol(beta)
  df1 = data.frame(cbind(lambda, beta))
  names(df1) = c("lambda", paste("beta", 1:p, sep= ""))
  df1.m = melt(df1, "lambda")
  names(df1.m)  = c("lambda", "Coefficients", "Estimated.Coefficients")
  Estimated.Coefficients = df1.m[,"Estimated.Coefficients"]
  Coefficients = df1.m[,"Coefficients"]
  g1 = ggplot(data=df1.m, aes(x = lambda, y= Estimated.Coefficients, colour= Coefficients)) + geom_line()
  g1
}

plot.orderedLasso.cv <- function(x, ...) {
  par(mar = c(5, 5, 5, 1))
  yrang=range(c(x$cv.err-x$cv.se,x$cv.err+x$cv.se))
  plot(log(x$lamlist), x$cv.err, xlab="log(lambda)",
       ylab = "Cross-validation Error", type="n",ylim=yrang)
  axis(3, at = log(x$lamlist), labels = paste(x$nonzeros), srt = 90, adj = 0)
  mtext("Number of features", 3, 4, cex = 1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  error.bars(log(x$lamlist), x$cv.err - x$cv.se, x$cv.err + x$cv.se, width = 0.01, col = "darkgrey")
  points(log(x$lamlist), x$cv.err, col=2, pch=19)
  abline(v=log(x$lamhat), lty=3)
  abline(v=log(x$lamhat.1se), lty=3)
}


predict.orderedLasso = function(object, newdata, ...){
    n <- nrow(newdata)
    ordered = object$strongly.ordered
    if (class(object$beta) == "numeric"){
      if(is.null(object$b0)){
        b0 = 0
        if(ordered) b0.ordered = 0
      }else{
        b0 = object$b0
        if(ordered) b0.ordered = object$b0.ordered
      }
      yhat = as.vector(newdata %*% object$beta + b0)
      if(ordered) yhat.ordered = as.vector(newdata %*% object$beta.ordered + b0.ordered)
    }else{
      if(is.null(object$b0)){
        b0  = rep(0, length = ncol(object$beta))
        if(ordered) b0.ordered  = rep(0, length = ncol(object$beta))
      }else{
        b0 = object$b0 
        if(ordered) b0.ordered  =object$b0.ordered
      }
      nlam = ncol(object$beta)
      yhat = matrix(NA, n, nlam)     
      if(ordered) yhat.ordered = matrix(NA, n, nlam)
      for (i in seq(nlam)){
        yhat[, i] = as.vector(newdata %*% object$beta[, i] + b0[i])
        if(ordered) yhat.ordered[ , i] = as.vector(newdata %*% object$beta.ordered[ ,i] + b0.ordered[i])
      }
    }
  if(!ordered) yhat.ordered = NULL
  return(list(yhat = yhat, yhat.ordered = yhat.ordered))
}

predict.orderedLasso.path = function(object, newdata, s= NULL, exact = FALSE, ...){
  ordered = object$strongly.ordered
  if(is.null(object$b0)){
    b0 = rep(0, length = ncol(object$beta))
    if(ordered) b0.ordered  = rep(0, length = ncol(object$beta))
  }else{
    b0 = object$b0
    if(ordered) b0.ordered = object$b0.ordered
  }
  if (!is.null(s))  s= sort(s ,decreasing = TRUE)
  if(missing(newdata)){
    stop("Missing data matrix")
  }
  lambda=object$lamlist
  if (exact && (!is.null(s))){
    which = match(s,lambda,FALSE)
    if(!all(which>0)){
      lambda_new = sort(s, decreasing = TRUE)
      object = update(object,lamlist=lambda_new)
      return (predict.orderedLasso(object = object, newdata = newdata, ...))
    }else{
      yhat = matrix(0, nrow = nrow(newdata), ncol = length(s))
      if(ordered) yhat.ordered = matrix(0, nrow = nrow(newdata), ncol = length(s))
   
      nbeta = rbind2((t(as.matrix(b0[which]))),object$beta[,which])
      if(ordered) nbeta.ordered = rbind((t(as.matrix(b0.ordered[which]))), object$beta.ordered[,which])
      
      for (i in (1 : ncol(nbeta))){
        temp = cbind(1, newdata)
        yhat[, i] =  temp %*% nbeta[,i]
        if(ordered) yhat.ordered[,i] = temp %*% nbeta.ordered[,i]
      }
      if(!ordered) yhat.ordered = NULL
      return(list(yhat = yhat, yhat.ordered = yhat.ordered))
    }
  }
  if(!is.null(s)){
    which = match(s,lambda,FALSE) 
    if(all(which>0)){
      yhat = matrix(0, nrow = nrow(newdata), ncol = length(s))
      if(ordered) yhat.ordered = matrix(0, nrow = nrow(newdata), ncol = length(s))
      b0 = t(as.matrix(b0[which]))
      if(ordered) b0.ordered = t(as.matrix(b0.ordered[which])) 
      nbeta = rbind2(b0,object$beta[,which])
      if(ordered) nbeta.ordered = rbind(b0.ordered, object$beta.ordered[,which])
      for (i in (1 : ncol(nbeta))){
        temp = cbind(1, newdata)
        yhat[, i] =  temp %*% nbeta[,i]
        if(ordered) yhat.ordered[,i] = temp %*% nbeta.ordered[,i]
      }
      if(!ordered) yhat.ordered = NULL
      return(list(yhat = yhat, yhat.ordered = yhat.ordered))
    }else{  
      b0 = t(as.matrix(b0))
      if(ordered) b0.ordered = t(as.matrix(b0.ordered))
      nbeta = rbind2(b0,object$beta)
      if(ordered) nbeta.ordered = rbind2(b0.ordered, object$beta.ordered)
      lambda = object$lamlist
      lamlist_interp = lambda.interp(lambda,s)
      yhat = matrix(0, nrow = nrow(newdata), ncol = length(s))
      if(ordered) yhat.ordered = matrix(0, nrow = nrow(newdata), ncol = length(s))
      nbeta = as.matrix(nbeta[,lamlist_interp$left,drop=FALSE]%*%Diagonal(x=lamlist_interp$frac) + 
                        nbeta[,lamlist_interp$right,drop=FALSE]%*%Diagonal(x=1-lamlist_interp$frac))    
      if(ordered) nbeta.ordered = as.matrix(nbeta.ordered[,lamlist_interp$left,drop=FALSE]%*%Diagonal(x=lamlist_interp$frac) + 
                        nbeta.ordered[,lamlist_interp$right,drop=FALSE]%*%Diagonal(x=1-lamlist_interp$frac))    
      for (i in (1 : ncol(nbeta))){
        temp = cbind(1, newdata)
        yhat[, i] =  temp %*% nbeta[,i]
        if(ordered) yhat.ordered[, i] =  temp %*% nbeta.ordered[,i]
      }
      if(!ordered) yhat.ordered = NULL
      return(list(yhat = yhat, yhat.ordered = yhat.ordered))
    }
  }
  predict.orderedLasso(object = object, newdata = newdata, ...)
}


timeLagLasso = function(x, y, lambda, maxlag, intercept = TRUE, standardize = TRUE, beta_pos = NULL, 
                        beta_neg = NULL, b0 = NULL, strongly.ordered = TRUE, maxiter = 500, inneriter = 100, iter.gg = 100, 
                        method = c("Solve.QP", "GG"), trace = FALSE, epsilon = 1e-5){
  ####error check
  stopifnot(nrow(x) == length(y), lambda >= 0)
  stopifnot(class(lambda) == "numeric")
  stopifnot(is.finite(x), is.finite(y), is.finite(lambda))
  if(missing(method)) method = "Solve.QP"

  this.call = match.call()

  if (standardize)  {
    stdeviation_x =  apply(x, 2, sd)
    stdeviation_inverse_scaled = rep(1/stdeviation_x, each = maxlag)
  }

  p = ncol(x)
  x = scale(x, center = FALSE, scale = standardize)

  if (is.null(beta_pos)) beta_pos = rep(0, (maxlag * ncol(x)))
  if (is.null(beta_neg)) beta_neg = rep(0, (maxlag * ncol(x)))
  cat("build the matrix for timeLagLasso, ", maxlag + 1, "observations in y are deleted \n")
  N = nrow(x)
  x = time_lag_matrix(x, maxlag)
  y = y[1: (N - maxlag - 1)]  
  est = timeLagLassoEstOrdered(x = x, y = y, lambda = lambda, strongly.ordered = strongly.ordered, 
                          maxlag = maxlag, intercept = intercept, b0 = b0,  beta_pos = beta_pos, 
                        beta_neg = beta_neg, stdeviation_inverse_scaled = stdeviation_inverse_scaled, standardize = standardize,
                        method = method, maxiter = maxiter, inneriter = inneriter, iter.gg = iter.gg, trace = trace, epsilon = epsilon)
  beta = est$beta
  beta_pos=  est$bp
  beta_neg = est$bn
  if(intercept) b0 = est$b0
  fitted = est$fitted
  err = est$err
 
  if(strongly.ordered){
    beta.ordered = est$beta.ordered
    if(intercept) b0.ordered =est$b0.ordered
    else b0.ordered = NULL
    err.ordered = est$err.ordered
    fitted.ordered = est$fitted.ordered
  }else{
    beta.ordered = b0.ordered = err.ordered = fitted.ordered = NULL
  }

  type = "gaussian"
  out = list(bp = beta_pos, bn = beta_neg, beta = beta, b0 = b0, 
             b0.ordered = b0.ordered, beta.ordered = beta.ordered, 
             maxlag = maxlag, p = p, strongly.ordered = strongly.ordered,
             fitted = fitted, fitted.ordered = fitted.ordered, err = err, err.ordered = err.ordered, 
             type = type,  call = this.call, method = method)
  class(out) = "timeLagLasso"
  out
}

timeLagLasso.path =function(x, y, lamlist = NULL, minlam = NULL, maxlam = NULL, nlam = 10, flmin = 1e-2, strongly.ordered = TRUE, 
                            flmax = 1, maxlag, intercept= TRUE, standardize = TRUE,  method = c("Solve.QP", "GG"),
                            maxiter = 500, inneriter = 100, iter.gg = 100, trace=FALSE, epsilon = 1e-5){  
  stopifnot(nrow(x) == length(y))
  stopifnot(is.finite(x), is.finite(y))
  this.call = match.call()
  if(missing(method)) method = "Solve.QP"
  if (is.null(maxlam)) {
    if (!is.null(minlam)) stop("Cannot have maxlam = NULL if minlam is non-null.")
    y_temp = y - mean(y)
    x_temp = scale(x, TRUE, TRUE)
    maxlam =  max(abs(crossprod(x_temp,y_temp)))
    minlam <- flmin * maxlam
  }

  if (is.null(minlam)) minlam <- maxlam * flmin
  if (is.null(lamlist))  lamlist <- exp(seq(log(maxlam), log(minlam),length=nlam))
  nlam <- length(lamlist)
  p = ncol(x)
  
  if (standardize){
    stdeviation_x =  apply(x, 2, sd)
    stdeviation_inverse_scaled = rep(1/stdeviation_x, each = maxlag)
  }
  ###
  x = scale(x, FALSE, scale = standardize) 
  if(trace) cat("build the matrix for timeLagLasso.path, ", maxlag + 1, "observations in y are deleted","\n")
  N = nrow(x)
  x = time_lag_matrix(x, maxlag) 
  y = y[1: (N - maxlag - 1)]
  ###
  beta_pos_lambda = matrix(0, ncol = nlam, nrow = (maxlag * p))
  beta_neg_lambda = matrix(0, ncol = nlam, nrow = (maxlag * p))
  beta = matrix(0, ncol = nlam, nrow = (maxlag * p)) 
  if(strongly.ordered) beta.ordered = matrix(0, ncol = nlam, nrow = (maxlag * p))
  bpinit = bninit = boinit = rep(0, maxlag*p)

  if(intercept){
    b0 = rep(NA, length = nlam)
    if(strongly.ordered)  b0.ordered = rep(NA, length = nlam)
    b0init = 0
  }else{
    b0 = b0init = NULL
    if(strongly.ordered) b0.ordered = NULL
  }
  fitted  = matrix(NA, ncol = nlam, nrow = nrow(x))
  err = rep(NA, length = nlam)
  
  if(strongly.ordered){
      fitted.ordered = matrix(NA, ncol = nlam, nrow = nrow(x))
      err.ordered = rep(NA, length = nlam)
  }

  for (i in 1: nlam){  
    if (trace) cat(c("lambda=",lamlist[i]),fill=T)
    if(i>1){
      bpinit = beta_pos_lambda[, i-1]
      bninit = beta_neg_lambda[, i-1]
      b0init = b0[i-1]
    }
    est = timeLagLassoEstOrdered(x = x, y = y, lambda = lamlist[i], maxlag = maxlag, intercept = intercept, b0 = b0init, 
                                 beta_pos=bpinit, beta_neg=bninit, strongly.ordered = strongly.ordered,
                                 stdeviation_inverse_scaled  = stdeviation_inverse_scaled, standardize = standardize, 
                                 method = method, maxiter = maxiter, inneriter = inneriter, trace= FALSE, epsilon = epsilon) 
    beta_pos_lambda[, i] = est$bp
    beta_neg_lambda[, i] = est$bn
    beta[,i] = est$beta
    if(intercept){
      b0[i] = est$b0 
    }
    fitted[, i] = est$fitted
    err[i] = est$err

    if(strongly.ordered){
       if(intercept) b0.ordered[i] = est$b0.ordered
       beta.ordered[,i] = est$beta.ordered
       fitted.ordered[, i ] = est$fitted.ordered  
       err.ordered[i] = est$err.ordered
    }else{
       beta.ordered = fitted.ordered = err.ordered = b0.ordered = NULL
    }   
  }
  out = list(bp = beta_pos_lambda, bn = beta_neg_lambda, beta = beta, b0 = b0, 
             beta.ordered = beta.ordered, 
             fitted.ordered = fitted.ordered, err.ordered = err.ordered,  b0.ordered = b0.ordered,  
             strongly.ordered = strongly.ordered, lamlist = lamlist, maxlag = maxlag, p = p,  
             fitted = fitted, err = err, call = this.call, method = method)
  class(out) = "timeLagLasso.path"
  out
}

timeLagLasso.cv = function(x, y, maxlag, lamlist = NULL, minlam = NULL, maxlam = NULL, nlam = 10, flmin = 1e-2, flmax = 1, 
                            intercept = TRUE, standardize = TRUE, method = c("Solve.QP", "GG"), strongly.ordered = TRUE, 
                            nfolds=10, folds = NULL,  maxiter = 500, inneriter = 100, iter.gg = 100, trace = FALSE, epsilon = 10e-5){
  # x: a matrix of predictors, where the rows are the samples and the columns are the predictors
  # y: a vector of observations, where length(y) equals nrow(x)
  # Split into nfolds, and for given folds, use the rest of data to get beta, beta_pos, beta_neg from  timeLagLassoPathEst 
  # and predict the current fold fitted values. 
  # Description for timeLagLassoPathEst 
  # Calculate the same thing with timeLagLasso.path except for taking well-build matrix x and y. 
  stopifnot(nrow(x) == length(y))
  stopifnot(is.finite(x), is.finite(y))
  if(missing(method)) method = "Solve.QP"
  p = ncol(x)
  errfun= function(yhat, y){(yhat - y)^2}
  this.call = match.call()

  if (is.null(maxlam)) {
    if (!is.null(minlam)) stop("Cannot have maxlam = NULL if minlam is non-null.")
    maxlam <- max(abs(crossprod(x,y)))
    minlam <- maxlam * flmin
    maxlam = maxlam * flmax
  }
  if (is.null(minlam)) minlam <- maxlam * flmin
  if (is.null(lamlist))  lamlist <- exp(seq(log(maxlam),log(minlam),length=nlam))
  
  nlam <- length(lamlist)
  fit = timeLagLasso.path(x, y, lamlist  = lamlist, minlam = minlam, maxlam = maxlam, nlam = nlam, flmin = flmin, 
                          flmax = flmax, maxlag = maxlag, intercept = intercept, strongly.ordered = strongly.ordered, 
                          standardize = standardize, maxiter  = maxiter, inneriter = inneriter, 
                          method = method, trace = FALSE, epsilon = epsilon)
  nonzeros = colSums(fit$beta !=0)
  if(strongly.ordered) nonzeros.ordered  = colSums(fit$beta.ordered != 0)
  else nonzeros.ordered = NULL

  if (standardize)  {
    stdeviation_x =  apply(x, 2, sd)
    stdeviation_inverse_scaled = rep(1/stdeviation_x, each = maxlag)
  }
  
  x = scale(x, FALSE, scale = standardize) 
  cat("build the matrix for timeLagLasso.cv, ", maxlag + 1, "observations in y are deleted\n")
  N = nrow(x)
  x = time_lag_matrix(x, maxlag)
  y = y[1: (N - maxlag - 1)]
  
  n <- length(y)
  if(is.null(folds)) {
    folds <- split(sample(1:n), rep(1:nfolds, length = n))
  }else {
    stopifnot(class(folds)=="list")
    nfolds <- length(folds)
  } 

  n.lamlist <- length(lamlist)   

  err2 = matrix(NA,nrow=nfolds,ncol=length(lamlist))
  if(strongly.ordered) err2.ordered = matrix(NA, nrow = nfolds, ncol = length(lamlist))

  for (ii in 1: nfolds){
    if(trace) cat("Fold", ii, ":\n")
    a = timeLagLassoPathEst(x = x[-folds[[ii]],], y=y[-folds[[ii]]], lamlist = lamlist, maxlag = maxlag, p = p, 
                            strongly.ordered = strongly.ordered, 
                            intercept = intercept, standardize = standardize, stdeviation_inverse_scaled = stdeviation_inverse_scaled, 
                            method = method, maxiter = maxiter, inneriter = inneriter, trace = trace)
    predict.model = predict.timeLagLassoPathEst(a, newdata=x[folds[[ii]],])
    yhatt = predict.model$yhat
    if(strongly.ordered) yhatt.ordered = predict.model$yhat.ordered
    temp = matrix(y[folds[[ii]]],nrow=length(folds[[ii]]),ncol=n.lamlist)
    err2[ii, ] = colMeans(errfun(yhatt,temp))
    if(strongly.ordered) err2.ordered[ii, ] = colMeans(errfun(yhatt.ordered, temp))
    if(trace) cat("\n")
  }
  errm=colMeans(err2)
  errse=sqrt(apply(err2,2,var)/nfolds)
  o=which.min(errm)
  lamhat=lamlist[o]
  oo=errm>= errm[o]+errse[o]
  lamhat.1se = sort(lamlist[oo & lamlist>=lamhat])[1]

  if(strongly.ordered){
    errm.ordered =colMeans(err2.ordered)
    errse.ordered=sqrt(apply(err2.ordered,2,var)/nfolds)
    o.ordered=which.min(errm.ordered)
    lamhat.ordered = lamlist[o.ordered]
    oo.ordered = errm.ordered >= errm.ordered[o.ordered] + errse.ordered[o.ordered]
    lamhat.ordered.1se = sort(lamlist[oo.ordered & lamlist>=lamhat.ordered])[1]
  }else{
    errm.ordered = errse.ordered = lamhat.ordered = lamhat.ordered.1se = NULL
  }
  obj <- list(lamlist = lamlist, cv.err = errm, cv.se=errse, lamhat=lamhat, lamhat.1se = lamhat.1se, 
              method = method,  strongly.ordered = strongly.ordered,
              nonzeros = nonzeros, nonzeros.ordered = nonzeros.ordered, 
              lamhat.ordered = lamhat.ordered, lamhat.ordered.1se = lamhat.ordered.1se, 
              cv.err.ordered = errm.ordered, cv.se.ordered = errse.ordered, call = this.call)
  class(obj) <- "timeLagLasso.cv"
  obj
}

#' plot coefficients from a "timeLagLasso.path" object
#' 
#' Produces a coefficient profile plot of the coefficient paths for a fitted "timeLagLasso.path" object.
#' @param x fitted "timeLagLasso.path" model
#' @param ...  Other graphical parameters to plot
#' @export

plot.timeLagLasso.path = function(x, ...){
  lambda = x$lamlist
  beta = t(x$beta)
  p = ncol(beta)
  df1 = data.frame(cbind(lambda = lambda, beta))
  names(df1) = c("lambda", paste("beta", 1:p, sep= ""))
  df1.m = melt(df1, "lambda")
  names(df1.m)  = c("lambda", "Coefficients", "Estimated.Coefficients")
  Coefficients = df1.m[,"Coefficients"]
  Estimated.Coefficients = df1.m[, "Estimated.Coefficients"]
  g1 = ggplot(df1.m, aes(x = lambda, y = Estimated.Coefficients, colour = Coefficients)) + geom_line()
  g1  
}

plot.timeLagLasso.cv <- function(x, ...) {
  par(mar = c(5, 5, 5, 1))
  yrang=range(c(x$cv.err-x$cv.se,x$cv.err+x$cv.se))
  plot(log(x$lamlist), x$cv.err, xlab="log(lambda)",
       ylab = "Cross-validation Error", type="n",ylim=yrang)
  axis(3, at = log(x$lamlist), labels = paste(x$nonzeros), srt = 90, adj = 0)
  mtext("Number of features", 3, 4, cex = 1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  error.bars(log(x$lamlist), x$cv.err - x$cv.se, x$cv.err + x$cv.se, width = 0.01, col = "darkgrey")
  points(log(x$lamlist), x$cv.err, col=2, pch=19)
  abline(v=log(x$lamhat), lty=3)
  abline(v=log(x$lamhat.1se), lty=3)
  
}

predict.timeLagLasso = function(object, newdata = NULL, ...){
  ordered = object$strongly.ordered
  if(is.null(newdata)){
    yhat = object$fitted
    if(ordered) yhat.ordered = object$fitted.ordered
  }else{
      if (class(object$beta) == "numeric"){
      if(is.null(object$b0)){
        b0 = 0
        if(ordered) b0.ordered = 0
      }else{
        b0 = object$b0
        if(ordered) b0.ordered = object$b0.ordered
      }
      yhat = as.vector(newdata %*% object$beta + b0)
      if(ordered) yhat.ordered =as.vector(newdata %*% object$beta.ordered + b0.ordered)
    }
    else{
      if(is.null(object$b0)){
        b0 = rep(0, length = length(object$beta))
        if(ordered) b0.ordered = rep(0, length = length(object$beta))
      }else{
        b0 = object$b0
        if(ordered) b0.ordered = object$b0.ordered
      }
      newdata = time_lag_matrix(newdata, object$maxlag)
      n <- nrow(newdata)
      nlam = length(object$lamlist)
      yhat = matrix(NA, nrow = n, ncol = nlam)
      if(ordered) yhat.ordered = matrix(NA, nrow = n, ncol = nlam)
      for (i in seq(nlam)){
       yhat[, i] = as.vector(newdata %*% object$beta[, i] + b0[i])
       if(ordered) yhat.ordered[, i] =as.vector(newdata %*% object$beta.ordered[ , i] + b0.ordered[i]) 
      }
    }
  }
  if(!ordered) yhat.ordered = NULL
  return (list(yhat = yhat, yhat.ordered = yhat.ordered))
}

predict.timeLagLasso.path = function(object, newdata = NULL, s = NULL, exact = FALSE, ...){
  ordered = object$strongly.ordered
  if (is.null(object$b0)){
    b0 = rep(0, length = ncol(object$beta))
    if(ordered) b0.ordered = rep(0, length = ncol(object$beta))
  }else{
    b0 = object$b0
    if(ordered) b0.ordered = object$b0.ordered
  }
  if (!is.null(s))  s= sort(s, decreasing = TRUE)
  if(missing(newdata)){
    stop("Missing data matrix")
  }

  lambda = object$lamlist
  if (exact && (!is.null(s))){
    which = match(s,lambda,FALSE)
    if(!all(which>0)){
      lambda_new = sort(s, decreasing = TRUE)
      object=update(object,lamlist=lambda_new)
      return (predict.timeLagLasso(object = object, newdata = newdata, ...))
    }else{
      newdata = time_lag_matrix(newdata, maxlag = object$maxlag)
      nfit = matrix(0, nrow = nrow(newdata), ncol = length(s))
      b0 = t(as.matrix(b0[which]))
      nbeta = rbind2(b0,object$beta[,which])
      temp = cbind(1, newdata)
      for (i in (1 : ncol(nbeta))){
        nfit[, i] =  temp %*% nbeta[,i]
      }
      if(ordered){
        nfit.ordered = matrix(0, nrow = nrow(newdata), ncol = length(s))
        b0.ordered = t(as.matrix(b0.ordered[which]))
        nbeta.ordered = rbind2(b0.ordered, object$beta.ordered[,which])
        for(i in 1:ncol(nbeta.ordered)){
          nfit.ordered = temp %*% nbeta.ordered[,i]  
        }
      }
      if(!ordered) nfit.ordered = NULL
      return(list(yhat = nfit, yhat.ordered = nfit.ordered))
    }
  }
  if(!is.null(s)){
     which = match(s,lambda,FALSE) 
     if(all(which>0)){
        newdata = time_lag_matrix(newdata, maxlag = object$maxlag)
        nfit = matrix(0, nrow = nrow(newdata), ncol = length(s))
        b0 = t(as.matrix(b0[which]))
        nbeta = rbind2(b0,object$beta[,which])
        temp = cbind(1, newdata)
        for (i in (1 : ncol(nbeta))){
           nfit[, i] =  temp %*% nbeta[,i]
        }
        
        if(ordered){
          nfit.ordered = matrix(0, nrow = nrow(newdata), ncol = length(s))
          b0.ordered = t(as.matrix(b0.ordered[which]))
          nbeta.ordered = rbind2(b0.ordered, object$beta.ordered[,which])
          temp = cbind(1, newdata)
          for (i in (1:ncol(nbeta))){
            nfit.ordered[,i] = temp %*% nbeta.ordered[, i]
          }
        } 
        if(!ordered) nfit.ordered = NULL
        return(yhat = nfit, yhat.ordered = nfit.ordered)
     }else{
        newdata = time_lag_matrix(newdata, maxlag = object$maxlag)
        b0 = t(as.matrix(b0))
        nbeta = rbind2(b0,object$beta)
        lambda = object$lamlist
        lamlist_interp = lambda.interp(lambda,s)
        nfit = matrix(0, nrow = nrow(newdata), ncol = length(s))
        nbeta = as.matrix(nbeta[,lamlist_interp$left,drop=FALSE]%*%Diagonal(x=lamlist_interp$frac) + 
                        nbeta[,lamlist_interp$right,drop=FALSE]%*%Diagonal(x=1-lamlist_interp$frac))    
        temp = cbind(1, newdata)
      
        for (i in (1 : ncol(nbeta))){
          nfit[, i] =  temp %*% nbeta[,i]
        }
        
        if(ordered){
          b0.ordered = t(as.matrix(b0.ordered))
          nbeta.ordered = rbind2(b0.ordered, object$beta.ordered)
          nfit.ordered = matrix(0, nrow = nrow(newdata), ncol = length(s))
          nbeta.ordered = as.matrix(nbeta.ordered[, lamlist_interp$left, drop = FALSE] %*% Diagonal(x = lamlist_interp$frac)+
                                    nbeta.ordered[, lamlist_interp$right, drop = FALSE] %*% Diagonal(x = 1 - lamlist_interp$frac))
          for(i in (1:ncol(nbeta.ordered))){
            nfit.ordered[,i] = temp %*% nbeta.ordered[,i]
          }
        }
      if(!ordered) nfit.ordered = NULL
      return(list(yhat = nfit, yhat.ordered = nfit.ordered))
    }
  }
  predict.timeLagLasso(object = object, newdata = newdata,...)
}

print.timeLagLasso = function(x, ...){
  cat("Call:\n")
  dput(x$call)
  if(x$strongly.ordered){
    mat = cbind(x$bp, x$bn, x$beta, x$beta.ordered)
    dimnames(mat) = list(NULL, c("beta_pos", "beta_neg", "beta", "beta.ordered"))

  }else{
    mat = cbind(x$bp, x$bn, x$beta)
    dimnames(mat) = list(NULL, c("beta_pos", "beta_neg", "beta"))
  }
  cat("\n")
  print(mat, quote = FALSE)
  cat("\n")

}

print.timeLagLasso.path <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
  mat = cbind(round(x$lamlist,2), round(x$err,3))
  dimnames(mat) = list(NULL, c("Lambda", "fittedError"))
  cat("\n")
  print(mat, quote = FALSE)
  cat("\n")
}

print.timeLagLasso.cv <- function(x, ...) {
  cat("Call\n")
  dput(x$call)
  if (x$strongly.ordered){
    mat = cbind(x$lamlist, x$cv.err, x$cv.se, x$cv.err.ordered, x$cv.se.ordered)
    dimnames(mat) = list(NULL, c("Lambda", "Mean CV Error", "SE", "Mean CV Error Ordered", "SE Ordered"))
  }else{
    mat = cbind(x$lamlist, x$cv.err, x$cv.se) 
    dimnames(mat) = list(NULL, c("Lambda", "Mean CV Error", "SE"))
  }
  cat("\n")
  print(mat, quote = FALSE)
  cat("\n")
  if(x$strongly.ordered){ 
    cat(c("lamhat=",round(x$lamhat,5),"lamhat.1se=",round(x$lamhat.1se,5), "lamhat.ordered = ", round(x$lamhat.ordered,5), "lamhat.ordered.1se = ", round(x$lamhat.ordered.1se,5)),fill=TRUE)
  }else{
    cat(c("lamhat=",round(x$lamhat,5),"lamhat.1se=",round(x$lamhat.1se,5)),fill = TRUE)
  }
 }

predict.timeLagLassoPathEst= function(object, newdata = NULL,...){
  predict.timeLagLassoInternal(object = object, newdata = newdata,...)
}

predict.timeLagLassoInternal = function(object, newdata, ...){
  ordered = object$strongly.ordered
  
  if(is.null(newdata)){
    yhat = object$fitted
    if(ordered) yhat.ordered = object$fitted.ordered
  }else{
      if (class(object$beta) == "numeric"){
      if(is.null(object$b0)){
        b0 = 0
        if(ordered) b0.ordered = 0
      }else{
        b0 = object$b0
        if(ordered) b0.ordered = object$b0.ordered
      }
      yhat = as.vector(newdata %*% object$beta + b0)
      if(ordered) yhat.ordered =as.vector(newdata %*% object$beta.ordered + b0.ordered)
    }
    else{
      if(is.null(object$b0)){
        b0 = rep(0, length = length(object$beta))
        if(ordered) b0.ordered = rep(0, length = length(object$beta))
      }else{
        b0 = object$b0
        if(ordered) b0.ordered = object$b0.ordered
      }
      n <- nrow(newdata)
      nlam = length(object$lamlist)
      yhat = matrix(NA, nrow = n, ncol = nlam)
      if(ordered) yhat.ordered = matrix(NA, nrow = n, ncol = nlam)
      for (i in seq(nlam)){
       yhat[, i] = as.vector(newdata %*% object$beta[, i] + b0[i])
       if(ordered) yhat.ordered[, i] =as.vector(newdata %*% object$beta.ordered[ , i] + b0.ordered[i]) 
      }
    }
  }
  if(!ordered) yhat.ordered = NULL
  return (list(yhat = yhat, yhat.ordered = yhat.ordered))
}

fastgg = function(x,y, sig, beta_pos, beta_neg, lam, t = 1, gam = 0.8, epsilon = 1e-5, trace = FALSE, inneriter.gg = 100){
  if(sig == 1) b = beta_pos
  if(sig == -1) b = beta_neg
  b_old_old = b
  b_old = b
  k = 1
  ii = 0
  while (TRUE & ii < inneriter.gg){
    ii = ii + 1
    alpha = b_old + (k-2)/(k+1) * (b_old - b_old_old)
    leastsquareconst =as.numeric(y - (x %*% alpha))
    const = -crossprod(x, leastsquareconst)
    constobjalpha = (1/2) * leastsquareconst %*% leastsquareconst  
    const3 = crossprod(const, alpha)
    junk = alpha - t * const
    b = prox(junk, t * lam)
    residualb = as.numeric((y - x %*% b))

    g0 = (1/2) * (residualb %*% residualb)
    g1 = constobjalpha + crossprod(const,b) - const3  + (1/(2 * t)) * ((b-alpha)%*%(b-alpha))
    diff = abs(g0 - g1)
    
    while (g0 > g1){
      t = gam * t
      junk = alpha - t * const
      b = prox(junk, t * lam)    
      h = as.numeric((y - x %*% b))
      g0 = (1/2) * (h %*% h)
      g1 = constobjalpha + crossprod(const,b) - const3 + (1/(2 * t)) * ((b-alpha)%*%(b-alpha))
    }
    if (sum(abs(b - b_old))< epsilon) break 
    b_old_old = b_old
    b_old = b
    k = k + 1  
  }
  if(ii == inneriter.gg) cat("generalized gradient failed to converge")
  if(sig ==  1) beta_pos = b
  if(sig == -1) beta_neg = b
  
  return(list(beta_pos=beta_pos,beta_neg=beta_neg))
}
ObjMinFun=function(x,y,bp,bn,lambda){
  # objective function
  yhat= x %*%(bp-bn) 
  (1/2)*sum( (y-yhat)^2)+lambda*sum(bp+bn)
}


ObjMinSignFun = function(x, y, b0.ordered, beta.ordered, lambda, sign){
  yhat = x %*% beta.ordered + b0.ordered
  (1/2) * sum((y - yhat)^2) + lambda * sum(sign * beta.ordered)
}

timeLagLassoEstOrdered = function(x, y, lambda, maxlag, intercept, b0, beta_pos, beta_neg, stdeviation_inverse_scaled, standardize,
                           method, strongly.ordered, maxiter = 500, inneriter = 100, iter.gg = 100, trace = TRUE, epsilon = 1e-5){
  # Description for timeLagLassoEst with "GG" algorithm (three method choices)
  # timeLagLassoEst takes a well-build matrix x_new, y, and minimize 1/2(y - x_new * (bp-bn))^2 + lambda*(bp-bn) st bp bn ge 0. 
  # p represents p predictors in the original matrix x. For given j = 1,...p, we update the corresponding beta_pos and beta_neg
  # and fix the rest of predictors constant by calling orderedLasso
  p = ncol(x)/maxlag
  beta_pos_new = beta_pos
  beta_neg_new = beta_neg 
  beta = beta_pos - beta_neg
  beta_new = beta
  ii=0
  go=TRUE

  if(intercept){
    mean_x = apply(x,2,mean)
    x = scale(x, TRUE, FALSE)
    mean_y = mean(y)
    y = y - mean_y
  }

  val =  ObjMinFun(x = x, y = y, bp = beta_pos_new, bn = beta_neg_new, lambda = lambda)      
  val_new = val
  if(trace)   cat(c("iteration", ii, "the objective value is ", val), fill = TRUE)
  
  while(go & ii< maxiter){  
    ii = ii + 1
    go = FALSE
    for (j in (0: (p-1))){      
      if (intercept) b0  = as.numeric(mean_y - mean_x %*% beta_new)
      else b0 = 0 
      subsetindex = (j * maxlag + 1) : (j * maxlag + maxlag)
      r_j =  y  - (x[, -subsetindex] %*% (beta_pos_new[-subsetindex] -  beta_neg_new[-subsetindex]))  
      beta_cal = orderedLasso(x[, subsetindex], r_j, lambda, intercept = FALSE, standardize = FALSE, 
                                      beta_pos = beta_pos[subsetindex], beta_neg = beta_neg[subsetindex], 
                                      method = method, trace = FALSE, niter = inneriter, iter.gg = iter.gg, strongly.ordered = FALSE)  
      beta_pos_new[subsetindex] = beta_cal$bp
      beta_neg_new[subsetindex] = beta_cal$bn
      beta_new[subsetindex] = beta_cal$beta  
    }
    val_new =  ObjMinFun(x = x, y = y, bp = beta_pos_new, bn = beta_neg_new, lambda = lambda)      
    if(trace) cat(c("iteration", ii, "the objective value is ", val), fill = TRUE)
    if (abs((val_new - val)/val) > epsilon) go = TRUE
    val = val_new
    beta_pos = beta_pos_new
    beta_neg= beta_neg_new
    beta = beta_new   
  }
  fitted = x %*% beta 
  err = mean((y - fitted)^2)
  ####strongly.ordered part#########################################################################################################
  b0.ordered  = b0
  temp = beta
  ordered.test = FALSE
  for(t in (0:(p-1))){
    index = (t * maxlag+1):(t*maxlag+maxlag)
    beta.temp = beta[index]
    beta.diff = beta.temp[1:(length(beta.temp)-1)]-beta.temp[2:(length(beta.temp))] 
    if(any(beta.diff < -1e-5)) ordered.test = TRUE
  }
  if(!intercept) b0 = NULL  
  if(strongly.ordered & ordered.test){  
    beta_ordered = beta_ordered_new = temp
    signvec = sign(temp)
    b0.ordered  = b0
    val_ordered = ObjMinSignFun(x = x, y = y, b0.ordered = 0, beta.ordered = beta_ordered, lambda = lambda, sign = signvec)
    jj = 1
    go = TRUE
    while(go & jj< maxiter){ 
      jj = jj + 1
      go = FALSE
      for (j in (0: (p-1))){  
       if (intercept) b0.ordered  = as.numeric(mean_y - mean_x %*% beta_ordered_new)
        else b0.ordered = 0 
        subsetindex = (j * maxlag + 1) : (j * maxlag + maxlag)
        r_j =  y   - (x[, -subsetindex] %*% (beta_ordered_new[-subsetindex]))
        beta_ordered_new[subsetindex] = ordLasSignPos(x = x[, subsetindex], y = r_j , lam = lambda, signvec = signvec[subsetindex])$bp
          
      }
      val_ordered_new =  ObjMinSignFun(x = x, y = y, b0.ordered = 0, beta.ordered = beta_ordered_new, lambda = lambda, sign  = signvec)     
      if(trace) cat(c("iteration", jj, "the second objective optimization value is ", val_ordered), fill = TRUE)
      if (abs((val_ordered_new - val_ordered)/val) > epsilon) go = TRUE
      val_ordered = val_ordered_new
      beta_ordered = beta_ordered_new
    }
  }
  if(strongly.ordered & (!ordered.test)) beta_ordered = beta
  if(!intercept) b0.ordered = NULL
  if(strongly.ordered) fitted.ordered = x %*% beta_ordered
  if(strongly.ordered) err.ordered = mean((y - fitted.ordered)^2)
  if (standardize){
    beta = stdeviation_inverse_scaled  * beta
    beta_pos =stdeviation_inverse_scaled * beta_pos
    beta_neg = stdeviation_inverse_scaled * beta_neg
    if(strongly.ordered) beta_ordered = stdeviation_inverse_scaled * beta_ordered
    else beta_ordered = NULL
  }
  if(!strongly.ordered) {
    beta_ordered = rep(NA, length(beta_pos))
    fitted.ordered = rep(NA, length(y))
    err.ordered = rep(NA, length(err))
  }

  return (list(bp = beta_pos, bn = beta_neg, beta = beta, b0 = b0, b0.ordered = b0.ordered, strongly.ordered = strongly.ordered,
               fitted = fitted, beta.ordered = beta_ordered, fitted.ordered = fitted.ordered, err = err, err.ordered = err.ordered))
}

timeLagLassoPathEst =function(x, y, lamlist, maxlag, p, intercept, standardize, method, stdeviation_inverse_scaled, maxiter,
                              inneriter, strongly.ordered, trace=FALSE){  
  beta_pos_lambda = matrix(0, ncol = length(lamlist), nrow = (maxlag * p))
  beta_neg_lambda = matrix(0, ncol = length(lamlist), nrow = (maxlag * p))
  beta = matrix(0, ncol = length(lamlist), nrow = (maxlag * p))
  if(strongly.ordered) beta.ordered = matrix(0, ncol = length(lamlist), nrow = (maxlag * p))
  bpinit = boinit = bninit = rep(0, maxlag*p)
  if(intercept){
    if(strongly.ordered) b0.ordered = rep(0, length = length(lamlist))
    b0 = rep(0, length = length(lamlist))
    b0init = 0
  }else{
    b0 = b0init = b0.ordered  = NULL
  }

  fitted  = matrix(NA, ncol = length(lamlist), nrow = nrow(x))
  if(strongly.ordered) fitted.ordered  = matrix(NA, ncol = length(lamlist), nrow = nrow(x))
  err = rep(NA, length = length(lamlist))
  if(strongly.ordered) err.ordered = rep(NA, length = length(lamlist))
  
  for (i in 1: length(lamlist)){
    if (trace) cat(c("lambda=",lamlist[i]),fill=T)
    if(i>1) {
      bpinit = beta_pos_lambda[, i-1]
      bninit = beta_neg_lambda[, i-1]
    }
    est = timeLagLassoEstOrdered(x = x, y = y, lambda = lamlist[i], intercept = intercept, 
                                 stdeviation_inverse_scaled = stdeviation_inverse_scaled, standardize = standardize, strongly.ordered = strongly.ordered, 
                                 maxlag = maxlag, b0 = b0init, beta_pos=bpinit, beta_neg=bninit, method = method, trace= FALSE, 
                                 maxiter = maxiter, inneriter = inneriter) 
    beta_pos_lambda[, i] = est$bp
    if(strongly.ordered) beta.ordered[, i] = est$beta.ordered
    beta_neg_lambda[, i] = est$bn
    beta[,i] = est$beta
    
    if(intercept) {
      b0[i] = est$b0
      if(strongly.ordered) b0.ordered[i] = est$b0.ordered
    }
  }
  if(!strongly.ordered) beta.ordered = b0.ordered = NULL
  out = list(bp = beta_pos_lambda, bn = beta_neg_lambda, beta = beta, beta.ordered = beta.ordered, b0 = b0, 
              b0.ordered = b0.ordered, lamlist = lamlist, maxlag = maxlag, p =  p, strongly.ordered = strongly.ordered)
  class(out) = "timeLagLassoPathEst"
  out
}

lambda.interp=function(lambda,s){
### lambda is the index sequence that is produced by the model
### s is the new vector at which evaluations are required.
### the value is a vector of left and right indices, and a vector of fractions.
### the new values are interpolated bewteen the two using the fraction
### Note: lambda decreases. you take:
### sfrac*left+(1-sfrac*right)
  if(length(lambda)==1){# degenerate case of only one lambda
    nums=length(s)
    left=rep(1,nums)
    right=left
    sfrac=rep(1,nums)
  }
  else{
      s[s > max(lambda)] = max(lambda)
      s[s < min(lambda)] = min(lambda)
      k=length(lambda)
      sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
      lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
      coord <- approx(lambda, seq(lambda), sfrac)$y
      left <- floor(coord)
      right <- ceiling(coord)
      sfrac=(sfrac-lambda[right])/(lambda[left] - lambda[right])
      sfrac[left==right]=1
    }
list(left=left,right=right,frac=sfrac)
}

print.orderedLasso = function(x, ...){
  cat("Call\n")
  dput(x$call)
  if(!is.null(x$beta.ordered)){
    mat = cbind(x$bp, x$bn, x$beta, x$beta.ordered )
    dimnames(mat) = list(NULL, c("beta_pos", "beta_neg", "beta", "beta.ordered"))
 
  }else{
    mat = cbind(x$bp, x$bn, x$beta)
    dimnames(mat) = list(NULL, c("beta_pos", "beta_neg", "beta"))
 
  }
  cat("\n")
  print(mat, quote = FALSE)
  cat("\n")
  
  if (!is.null(x$b0)){
    if(x$strongly.ordered){
      cat("beta_zero =", x$b0, "beta.ordered_zero = ", x$b0.ordered)
    }else{
      cat("beta_zero =", x$b0)
    }
    cat("\n")
  }
}

print.orderedLasso.path = function(x,...){
  cat("Call\n")
  dput(x$call)
  ordered = x$strongly.ordered
  if(ordered){
    mat = cbind(x$lamlist, x$err, x$err.ordered)
    dimnames(mat) = list(NULL, c("Lambda", "Error", "Error.ordered"))
  }else{
    mat = cbind(x$lamlist, x$err)
    dimnames(mat) = list(NULL, c("Lambda", "Error"))
  }
  cat("\n")
  print(mat, quote = FALSE)
  cat("\n")
}

error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
  invisible()
}

generalizedGradient = function(x, y, sig, beta_pos, beta_neg, b0 = NULL, lambda, 
                               tt = 1, trace=FALSE, adjustb=NULL, backtrack= NULL, epsilon = 1e-12){
  # generalized gradient algorithm for isotonic coef reg model
  # sig=1 means we're solving for beta_pos;  -1 means beta_neg
  # we get update beta and step from backtrack function and do this until convergence. (since the step tt is calculated from backtrack
  #this function does not take a lot of iterations and thus I removed the number of iterations limit)
  if(sig == 1) b = beta_pos
  if(sig == -1) b = beta_neg
  val = ObjMinFun(x,y, beta_pos, beta_neg, lambda)
  ii = 0
  while(TRUE & ii < 500){
    ii = ii + 1
    b_old = b
    update  = backtrackfun(b, lambda, x, y, tt)
    tt = update$tt
    b =  update$beta
    if (sum(abs(b - b_old))< epsilon) break 
  } 
  if(sig ==  1) beta_pos = b
  if(sig == -1) beta_neg = b
  return(list(beta_pos=beta_pos,beta_neg=beta_neg))
}

backtrackfun=function(bp,lam,x,y,tt=1,gam=.8,trace=F){
  bn = rep(0,length(bp))
  g0 = ObjMinFun(x,y,bp,bn,0)
  dg = -t(x)%*%(y-x%*%bp)
  go = T
  i = 0
  while(go && i < 200){
    i = i + 1
    go = F
    r = bp-tt*dg
    G = (1/tt)*(bp-prox(r, tt*lam))
    val = g0-tt* sum(dg* G) +(tt/2)*sum(G^2)
    beta1 = bp-tt*G
    g1 = ObjMinFun(x,y,beta1,bn,0) 
    if(g1 > val){
      go=T
      tt=gam*tt
    }
    if(trace){cat(c(tt,val,g1),fill=T)}
  }
  return(list(beta=beta1,tt=tt))
}

ObjMinFun=function(x,y,bp,bn,lambda){
  yhat= x %*%(bp-bn) 
  (1/2)*sum( (y-yhat)^2)+lambda*sum(bp+bn)
}

prox=function(y,lam){
  aa=pava(y-lam, decreasing=TRUE)
  bb=aa*(aa>0)
}
#' @title time_lag_matrix
#' @description Build the time lag matrix for the data matrix x.
#' @param x A matrix of predictors, where the rows are the samples and the columns are the predictors
#' @param maxlag maximum time lag variable.
#' @export
time_lag_matrix = function (x, maxlag){
  p = ncol(x)
  N = nrow(x) - maxlag - 1
  x_new_num_col = maxlag * p 
  x_new = matrix(0, nrow =  N, ncol = x_new_num_col) 
  for (i in (1: N)){
    x_temp = x[(i + 1):(i + maxlag), ]
    x_new[i, ] = as.vector(x_temp)
  }
  return (x_new)
} 

ordLasSignPos = function(x,y,lam, signvec){
  n = nrow(x)
  A=diag(signvec)
  p = ncol(x)
  A[col(A)==(row(A)+1)]=-signvec[2:length(signvec)]
  A=A[-p,,drop=F]
  AA=rbind(A,diag(signvec))
  b=rep(0,nrow(AA))
  Amat=t(AA)
  emin=-min(eigen(t(x)%*%x)$val)
  emin=max(.001,emin)
  dmat=t(x)%*%x+diag(rep(emin,p))
  dvec=as.vector((t(x)%*%y)- lam * signvec)
  bvec = rep(0, length = nrow(AA))
  a=solve.QP(dmat, dvec, Amat, bvec)
  bp=a$sol
  invisible()
  beta.temp = bp*signvec
  return(list(bp=bp))
}

ordLas2=function(x,y,lam){
  # compute the ordered lasso using solve.QP
  # min .5*sum( y-x*(bp-bn))^2 +lam*sum(bp+bn)  st bp, bn ge 0
  n=nrow(x)
  xx=cbind(x,-x)
  p=ncol(x)
  BIG=10e9
  A=diag(p)
  A[col(A)==(row(A)+1)]=-1
  A=A[-p,, drop = F]
  AA=matrix(0,nrow=2*(p-1),ncol=2*p)
  AA[1:(p-1),1:p]=A
  AA[p:(2*(p-1)),(p+1):(2*p)]=A
  AA=rbind(AA,diag(2*p))
  b=rep(0,nrow(AA)) 
  r=rep(BIG,nrow(AA))
  Amat=t(AA)
  # had to add a small constant, since solve.QP wants dmat to be PD
  emin=-min(eigen(t(xx)%*%xx)$val)
  emin=max(.001,emin)
  dmat=t(xx)%*%xx+diag(rep(emin,2*p))
  dvec=as.vector(t(xx)%*%y -  rep(lam,2*p))
  bvec = rep(0, length = nrow(AA))
  a=solve.QP(dmat, dvec, Amat, bvec)
  bp=a$sol[1:p]
  bn=a$sol[-(1:p)]
  return(list(b=bp-bn,bp=bp,bn=bn))
}