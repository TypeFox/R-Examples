predict.ibr <- function(object, newdata, interval= c("none", "confidence", "prediction"), ...) {
    if (!inherits(object, "ibr")) stop("must be an ibr object\n")
    interval <- match.arg(interval)
    if ((interval == "prediction")|(interval == "confidence")) {
      warning("Interval for prediction/confidence is not implemented yet\n")
    }
  if (missing(newdata) || is.null(newdata)) {
    Yres <- object$fitted
  } else {
      if (any(is.na(newdata))) stop("NA's in newdata\n")
      tt <- terms(object)
      Terms <- delete.response(tt)
      m <- model.frame(Terms, newdata, na.action = na.fail)
      newdata <- model.matrix(Terms, m)
      attributes(newdata) <- attributes(newdata)[c("dim","dimnames")]
      t <- terms(object)
      mf <- object$call
      m <- match(c("formula", "data", "subset"), names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval(mf, environment(t))
      x <- model.matrix(t, mf)
    if (object$parcall$scaled) {
        newdata  <- scale(newdata,center=object$parcall$mean,scale=object$parcall$sd)
        x <- scale(x,center=object$parcall$mean,scale=object$parcall$sd)
    }
    if (object$parcall$smoother=="k") {
      SSx <- kernelSx(kernelx=object$parcall$kernel,X=x,bx=object$bandwidth,newdata)
      Yres <- as.vector(SSx%*%object$beta)
    }
    if ((object$parcall$smoother=="ds")|(object$parcall$smoother=="tps")) {
      SSx <- dsSx(x,newdata,m=object$parcall$m,s=object$parcall$s)
      Yres <- as.vector(SSx$Sgu%*%object$beta$d+SSx$Qgu%*%object$beta$c)
    }    
    if ((object$parcall$smoother=="lrds")|(object$parcall$smoother=="lrtps")) {
        SSx <- PredictMat(object$parcall$smoothobject,data.frame(newdata))
        Yres <- as.vector(SSx%*%object$beta)
    }    
  }
  return(Yres)
}

