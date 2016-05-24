RMSEP <- function(object, ...) UseMethod("RMSEP")

RMSEP.default <- function(object, ...)
  {
    stop("No default method for \"RMSEP\"")
  }

RMSEP.mat <- function(object, k, weighted = FALSE, ...) {
  if(!inherits(object, "mat"))
    stop("'object' is not of class \"mat\".")
  if(missing(k))
    k <- getK(object)
  if(weighted)
    rmsep <- object$standard$rmsep[k]
  else
    rmsep <- object$weighted$rmsep[k]
  return(rmsep)
}

RMSEP.bootstrap.mat <- function(object, type = c("birks1990", "standard"), ...) {
  if(!inherits(object, "bootstrap.mat"))
    stop("'object' is not of class \"bootstrap.mat\".")
  if(missing(type))
    type <- "birks1990"
  type <- match.arg(type)
  if(type == "birks1990")
    rmsep <- object$bootstrap$rmsep[getK(object)]
  else
    rmsep <- sqrt(mean(object$bootstrap$residuals[, getK(object)]^2))
  return(rmsep)
}
