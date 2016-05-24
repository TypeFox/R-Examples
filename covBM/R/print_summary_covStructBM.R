#' print.summary.corStructBM
#'
#' @param x An object of class "summary.corStructBM", containing information
#' on fitted stochastic process component.
#' @param ... Additional arguments (not used for this method).
#' @export
print.summary.corStructBM <-
  function(x, ...)
{
  class(x) <- attr(x, "oClass")
  cat(paste("Stochastic process component: ", attr(x, "structName"), "\n", sep = ""))
  cat(paste(" Formula:", deparse(formula(x)),"\n"))
  cat(" Parameter estimate(s):\n")
  if(is.null(attr(x, "sigma"))) {stop("Brownian motion or IOU components must be fitted using 'lmeBM' or 'nlmeBM' functions")}
  sigma<-attr(x, "sigma")
  if(inherits(x, "covBM")){
    print(coef(x, unconstrained = FALSE)*(sigma^2))	###sigma is multiplied back in to give natural parameter
    }  else if (inherits(x, "covFracBM")) {
    coefNat<-coef(x, unconstrained = FALSE)
    coefNat[1]<-coefNat[1]*(sigma^2)
    print(coefNat)
    }  else if (inherits(x, "covIOU")) {
    coefNat<-coef(x, unconstrained = FALSE)
    coefNat[1]<-coefNat[1]*(sigma^2)
    print(coefNat)
    }	else {print(coef(x, unconstrained = FALSE))}
  invisible(x)
}

#' summary.covBM
#'
#' @param object An object of class "covBM", containing information
#' on fitted stochastic process component.
#' @param structName An optional character string defining the type of correlation
#' structure associated with object, as for \code{\link[nlme]{summary.corStruct}}. Defaults to class(object)[1].
#' @param ... Additional arguments (not used for this method).
#' @export
summary.covBM <-
  function(object, structName = class(object)[1], ...)
{
  attr(object, "structName") <- structName
  attr(object, "oClass") <- class(object)
  class(object) <- "summary.corStructBM"
  object
}

#' summary.covFracBM
#'
#' @param object An object of class "covFracBM", containing information
#' on fitted stochastic process component.
#' @param structName An optional character string defining the type of correlation
#' structure associated with object, as for \code{\link[nlme]{summary.corStruct}}. Defaults to class(object)[1].
#' @param ... Additional arguments (not used for this method).
#' @export
summary.covFracBM <-
  function(object, structName = class(object)[1], ...)
{
  attr(object, "structName") <- structName
  attr(object, "oClass") <- class(object)
  class(object) <- "summary.corStructBM"
  object
}

#' summary.covIOU
#'
#' @param object An object of class "covIOU", containing information
#' on fitted stochastic process component.
#' @param structName An optional character string defining the type of correlation
#' structure associated with object, as for \code{\link[nlme]{summary.corStruct}}. Defaults to class(object)[1].
#' @param ... Additional arguments (not used for this method).
#' @export
summary.covIOU <-
  function(object, structName = class(object)[1], ...)
{
  attr(object, "structName") <- structName
  attr(object, "oClass") <- class(object)
  class(object) <- "summary.corStructBM"
  object
}
