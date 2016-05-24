## functions for extracting and setting the number
## of analogues to be used
getK <- function(object, ...) UseMethod("getK")

getK.default <- function(object, ...) {
  stop("No default method for 'k'")
}

getK.mat <- function(object, weighted=FALSE, ...){
  ## check that this is a mat object
  if(class(object) != "mat")
    stop("'object' must be of class 'mat'")
  if(weighted){
    retval <- object$weighted$k
    attr(retval, "auto") <- object$weighted$auto
    attr(retval, "weighted") <- TRUE
  } else {
    retval <- object$standard$k
    attr(retval, "auto") <- object$standard$auto
    attr(retval, "weighted") <- FALSE
  }
  return(retval)
}

getK.bootstrap.mat <- function(object, which = c("bootstrap", "model"),
                               prediction = FALSE, ...) {
  if (!inherits(object, "bootstrap.mat"))
    stop("Use only with \"bootstrap.mat\" objects")
  if(missing(which))
    which <- "bootstrap"
  which <- match.arg(which)
  if(which == "bootstrap") {
    if(prediction)
      retval <- object$predictions$bootstrap$k
    else
      retval <- object$bootstrap$k
  } else {
    if(prediction)
      retval <- object$predictions$model$k
    else
      retval <- object$predictions$model$k
  }
  attr(retval, "auto") <- object$auto
  attr(retval, "weighted") <- object$weighted
  return(retval)
}

getK.predict.mat <- function(object,
                             which = c("model", "bootstrap"),
                             ...) {
    if(missing(which))
        which <- "model"
    which <- match.arg(which)
    if(which == "bootstrap" && is.null()) {
        which <- "model"
        warning()
    }
    if(which == "model") {
        retval <- object$predictions$model$k
    } else {
        retval <- object$predictions$bootstrap$k
    }
    attr(retval, "auto") <- object$auto
    attr(retval, "weighted") <- object$weighted
    return(retval)
}

"setK<-" <- function(object, weighted=FALSE, value) UseMethod("setK<-")

"setK<-.default" <- function(object, weighted=FALSE, value) {
  stop("no default replacement method for 'k'")
}

"setK<-.mat" <- function(object, weighted=FALSE, value) {
  ## check that this is a mat object
  if(class(object) != "mat")
    stop("'object' must be of class 'mat'")
  ## check that value is not NULL
  if(is.null(value))
    stop("attempt to set NULL number of analogues, 'k'")
  ## check that value is numeric
  ## need to correct this, is.numeric is not integer
  if(!is.numeric(value))
    stop("attempt to set non-integer number of analogues, 'k'")
  if(weighted) {
    object$weighted$k <- value
    object$weighted$auto <- FALSE
  } else {
    object$standard$k <- value
    object$standard$auto <- FALSE
  }
  object
}
