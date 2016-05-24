predict.npregress <- function(object, newdata, interval= c("none", "confidence", "prediction"), deriv=FALSE, ...) {
  interval <- match.arg(interval)
  if ((interval == "prediction")|(interval == "confidence")) {
    warning("Interval for prediction/confidence is not implemented yet\n")
  }
    x <- object$call$x
    y <- object$call$y
  if (missing(newdata) || is.null(newdata)) {
    if (!deriv) return(object$fitted) else newdata <- x
  } else {
    if (any(is.na(newdata))) stop("NA's in newdata\n")
    if (!is.numeric(newdata)&(is.data.frame(newdata))) {
      newdata <- newdata[,1]
      if (!is.numeric(newdata)) stop("first column of data-frame is not numeric\n")
    }
    if (is.matrix(newdata)) {
      newdata <- as.vector(newdata)
    }
    if (!is.numeric(newdata)) stop("newdata must be a numeric vector (or a data-frame with first column of numeric type)\n")
  }
  ## autre methode
  if (object$call$degree==0) {
    methode <- "reg"
    nom <- paste(methode,object$call$kernel,sep="")
    prov <- .C(nom,as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1))
    deriv <- FALSE
  }
  if (object$call$degree==1) {
    methode <- "regpol"
    nom <- paste(methode,object$call$kernel,sep="")
    prov <- .C(nom,as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1),double(length(newdata)))
  }
  if (!deriv) {
    Yres <- prov[[7]]
  } else {
    Yres <- list(yhat=prov[[7]],deriv=prov[[9]])  
  }
  if (object$call$degree>1) stop("Not implemented. Please consider using KernSmooth or another library for degree greater or equal to 2\n")
  return(Yres)
}

