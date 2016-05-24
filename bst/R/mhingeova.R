mhingeova <- function(xtr, ytr, xte=NULL, yte=NULL, cost = NULL, nu=0.1, learner=c("tree", "ls", "sm"), maxdepth=1, m1=200, twinboost=FALSE, m2=200){
  call <- match.call()
  n1 <- length(ytr)
  ncla <- length(unique(ytr))
  if(min(ytr)!=1 || max(ytr)!=ncla)
    stop("Response variable must be 1, 2, ...\n")
  if(!is.null(xte) && dim(xtr)[2] != dim(xte)[2])
    stop("The training data and test data have different dimensions\n")  
  if(is.null(cost))
    cost.mul <- rep(0.5, ncla)
  else {
    cost.mul <- cost
    if(any(cost <= 0) || any(cost >= 1))
      stop("The costs must in (0, 1)\n")
  }
  z <- matrix(-1, ncol=ncla, nrow=n1)
  if(!is.null(xte) && !is.null(yte))
    z.te <- matrix(-1, ncol=ncla, nrow=length(yte))
  fit.tr <- fit.te <- fpar <- ensemble <- vector("list", ncla)
  err.tr <- err.te <- risk.te <- xsel <- NULL
  for(k in 1:ncla){
    z[which(ytr==k), k] <- 1
    if(!is.null(xte) && !is.null(yte))
      z.te[which(yte==k), k] <- 1
    cost1 <- cost.mul[k]
    eval(parse(text = paste("tr",k, " <- bst(x = xtr, y = z[,k], cost=cost.mul[k], learner=learner, ctrl = bst_control(mstop=m1, nu=nu), control.tree=list(maxdepth=maxdepth))", sep=""))) 
    if(twinboost){
      eval(parse(text = paste("tr",k, " <- bst(x = xtr, y = z[,k],  cost=cost.mul[k], learner=learner, ctrl = bst_control(twinboost=TRUE, f.init=predict(tr", k, "), xselect.init = tr", k, "$xselect, mstop=m2, nu=nu), control.tree=list(maxdepth=maxdepth))", sep=""))) 
    } 
    fit.tr[[k]] <- eval(parse(text=paste("predict(tr",k,", type='all.res')", sep="")))           
    if(learner!="tree" || maxdepth==1) 
    fpar[[k]] <- eval(parse(text=paste("fpartial.bst(tr",k,")", sep="")))           
    xsel <- c(xsel, eval(parse(text=paste("tr",k,"$xselect", sep=""))))           
    ensemble[[k]] <- eval(parse(text=paste("tr",k,"$ensemble", sep="")))           
    if(k==1)
      mstop1 <- eval(parse(text=paste("tr",k,"$ctrl$mstop", sep="")))           
    if(!is.null(xte)){
      fit.te[[k]] <- eval(parse(text=paste("predict(tr",k,", newdata=xte, type='all.res')", sep=""))) 
      if(!is.null(yte)){
        risk.te <- cbind(risk.te, eval(parse(text=paste("predict(tr",k,", newdata=xte, newy=z.te[,k], type='loss')", sep=""))))           
      }
    }
  }
  err.tr <- unlist(lapply(1:mstop1, function(j){tmp <- NULL; for(i in 1:ncla) tmp <- cbind(tmp, fit.tr[[i]][,j]); 1/length(ytr) * sum(ytr != apply(tmp, 1, which.max)) }))
  if(!is.null(yte) && !is.null(xte)){
    err.te <- unlist(lapply(1:mstop1, function(j){tmp <- NULL; for(i in 1:ncla) tmp <- cbind(tmp, fit.te[[i]][,j]); 1/length(yte) * sum(yte != apply(tmp, 1, which.max)) }))
}
  RET=list(call=call, learner=learner, nu=nu, twinboost=twinboost, m1=m1, m2=m2, risk.te=risk.te, err.tr = err.tr, err.te = err.te, ensemble = ensemble, xsel=round(table(xsel)/ncla,2), fpar=fpar)
  class(RET) = "mhingeova"
  return(RET)
}

cv.mhingeova <- function(x, y, balance=FALSE, K=10, cost = NULL, nu=0.1, learner=c("tree", "ls", "sm"), maxdepth=1, m1=200, twinboost = FALSE, m2=200, trace=FALSE, plot.it = TRUE, se = TRUE, ...)
{
  call <- match.call()
  learner <- match.arg(learner)
  if(balance)
  all.folds <- balanced.folds(y, K)
  else all.folds <- cv.folds(length(y), K)
  fraction <- 1:m1
#  fraction <- seq(from = 1, to = m1, by=5)
  residmat <- matrix(0, length(fraction), K)
  for(i in seq(K)) {
    if(trace)
      cat("\n CV Fold", i, "\n\n")
    omit <- all.folds[[i]]
    fit <- mhingeova(xtr = x[ - omit,,drop=FALSE  ], ytr = y[ - omit], xte = x[ omit,,drop=FALSE ], yte=y[ omit ], cost = cost, nu=nu, learner = learner, maxdepth=maxdepth, m1=m1, twinboost=twinboost, m2=m2)
###cross validation misclassification error
    residmat[,i] <- fit$err.te
  if(trace && i==K)
  cat("End of cross-validation\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(residmat = residmat, fraction = fraction, cv = cv, cv.error = cv.error)
  if(plot.it) plotCVbst(object,se=se)
  invisible(object)
}

print.mhingeova <- function(x, ...) {

  cat("\n")
  cat("\t Multi-class HingeBoost Fitted with One-against-All\n")
  cat("\n")
  if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  cat("Base learner: ", x$learner, "\n")
  cat("Step size: ", x$nu, "\n")
  if(!x$cv1)
    cat("Number of boosting iterations: mstop =", x$m1, "\n")
  if(x$twinboost){
    cat("Twin boosting", "\n")
    if(!x$cv2)
      cat("Number of twin boosting iterations: mstop =", x$m2, "\n")
  }
  invisible(x)
}


