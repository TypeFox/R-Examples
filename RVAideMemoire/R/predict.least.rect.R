predict.least.rect <- function(object,newdata,...) {
  if (class(object)!="least.rect") {stop("incorrect model")}
  pred.fun <- function(x,lev=NULL) {
    if (is.null(lev)) {
	x*object$coefficients[2]+object$coefficients[1]
    } else {
	x*object$coefficients[lev,2]+object$coefficients[lev,1]
    }
  }
  predicted <- NULL
  if (class(newdata)=="list") {
    if (!all(colnames(object$model)[-1]%in%names(newdata))) {stop("incorrect variable names")}
    covar <- which(names(newdata)==colnames(object$model)[2])
    fact <- if (ncol(object$model)==3) {
	which(names(newdata)==colnames(object$model)[3])
    } else {NULL}
    if (length(newdata)>1 & !is.null(fact)) {
	if (length(newdata[[covar]])!=length(newdata[[fact]])) {stop("variable lengths differ")}
    }
    for (i in 1:length(newdata[[covar]])) {
	if (is.null(fact)) {
	  predicted <- c(predicted,pred.fun(newdata[[covar]][i]))
	} else {
	  predicted <- c(predicted,pred.fun(newdata[[covar]][i],as.character(newdata[[fact]][i])))
	}
    }
  } else if (class(newdata)=="data.frame") {
    if (!all(colnames(object$model)[-1]%in%colnames(newdata))) {stop("incorrect variable names")}
    covar <- which(colnames(newdata)==colnames(object$model)[2])
    fact <- if (ncol(object$model)==3) {
	which(colnames(newdata)==colnames(object$model)[3])
    } else {NULL}
    for (i in 1:nrow(newdata)) {
	if (is.null(fact)) {
	  predicted <- c(predicted,pred.fun(newdata[i,covar]))
	} else {
	  predicted <- c(predicted,pred.fun(newdata[i,covar],as.character(newdata[i,fact])))
	}
    }
  } else {stop("incorrect 'newdata' argument")}
  names(predicted) <- NULL
  return(predicted)
}
