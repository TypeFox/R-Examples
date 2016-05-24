#' \code{derivateMF} derivate membership function
#'
#' Derivate de membership of x with respect to i of MembershipFunction object 
#' heirs.
#' 
#' @param object MembershipFunction class heirs
#' @param x numeric of the MembershipFunction to be evaluated
#' @param i index of the ith parameter to partially derivate
#'
#' @return numeric with the value obtained from the ith derivative at x
#'
#' @include MembershipFunction-evaluateMF.R
#' @exportMethod derivateMF
#' @docType methods
#' @name derivateMF
#' @rdname derivateMF
#' @aliases derivateMF-methods
#' @seealso \code{\link{MembershipFunction-class}} and \code{\link{evaluateMF}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setGeneric(name="derivateMF", def=function(object, x, i){
  standardGeneric("derivateMF")
})
#'
#' @name derivateMF
#' @rdname derivateMF
#' @inheritParams derivateMF
#' @aliases derivateMF,MembershipFunction-method
setMethod(f="derivateMF", signature="MembershipFunction", 
  definition=function(object, x, i){
  stop("Only heir of MembershipFunction class can derivateMF!!")
})
#'
#' @name derivateMF
#' @rdname derivateMF
#' @aliases derivateMF,BellMF-method
#' @inheritParams derivateMF
#' @examples
#' #BellMF example I
#' #A bell membership function with default prototype (a=1, b=1,c=0)
#' #The membership of x in the bell, should be 1
#' #The derivate of the first parameter at x, should be 0
#' #The derivate of the first parameter at x, should be also 0
#' bell <- new(Class="BellMF")
#' bell
#' evaluateMF(object=bell, x=0)
#' derivateMF(object=bell, x=0, i=1)
#' derivateMF(object=bell, x=0, i="a")
#' #
#' #BellMF example II
#' #A bell membership function with parameters (a=4,b=1,c=-10)
#' #The membership of x in the bell, should be 0.137931
#' #The derivate of the first parameter at x, should be 0.05945303
#' #The derivate on "a" at x=0, should be 0.05945303
#' bell2 <- new(Class="BellMF",parameters=c(a=4,b=1,c=-10))
#' bell2
#' evaluateMF(object=bell2, x=0)
#' derivateMF(object=bell2, x=0, i=1)
#' derivateMF(object=bell2, x=0, i="a")
setMethod(f="derivateMF",signature = "BellMF", 
  definition = function(object,x,i){
  out <- switch(i,
    a=(2*(object@parameters["b"]^2)*(-object@parameters["c"] + x)^2*
      ((((-object@parameters["c"] + x)/
      object@parameters["a"])^2)^(-1 + object@parameters["b"]^2)))/
      (object@parameters["a"]^3*(1 + (((-object@parameters["c"] + x)/
      object@parameters["a"])^2)^(object@parameters["b"]^2))^2),
    b=-(2*object@parameters["b"]*((((-object@parameters["c"] + x)/
      object@parameters["a"])^2)^(object@parameters["b"]^2))*
      log(((-object@parameters["c"] + x)/object@parameters["a"])^2))/
      ((1+(((-object@parameters["c"]+x)/
      object@parameters["a"])^2)^(object@parameters["b"]^2))^2),
    c=(2*(object@parameters["b"]^2)*(-object@parameters["c"] + x)*
      (((-object@parameters["c"] + x)/object@parameters["a"])^2)^(-1 + 
      object@parameters["b"]^2))/((object@parameters["a"]^2)*(1 + 
      (((-object@parameters["c"] + x)/
      object@parameters["a"])^2)^(object@parameters["b"]^2))^2)
  )
  ##Check option
  if(is.null(out)){
    stop("Not correct index i for the partial derivate")
  }
  #For the undefined of ther derivative_c lim x+->c = 0
  if(is.nan(out) & x==object@parameters["c"] & (i==2 | i=="b")){
    out <- 0 
  }
  if(is.nan(out)){
    stop(paste("Partial derivate returned NaN for x=", x, " i=", i, 
      " parameters=", toString(object@parameters),sep="")) 
  }
  return(out)
})
#'
#' @name derivateMF
#' @rdname derivateMF
#' @aliases derivateMF,GaussianMF-method
#' @inheritParams derivateMF
#' @examples
#' #GaussianMF example I
#' #A Gaussian membership function with default prototype (mu=0, sigma=1)
#' #The membership of x in the gaussian, should be 1/sqrt(2*pi) = 0.3989423
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' gaussian <- new(Class="GaussianMF")
#' gaussian
#' evaluateMF(object=gaussian, x=0)
#' derivateMF(object=gaussian, x=0, i=1)
#' derivateMF(object=gaussian, x=0, i="mu")
#' #
#' #GaussianMF example II
#' #A Gaussian membership function with parameters (mu=0, sigma=1)
#' #The membership of x in the Gaussian, should be 1/sqrt(2*pi) = 0.3989423
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' gaussian2 <- new(Class="GaussianMF",parameters=c(mu=0,sigma=1))
#' gaussian2
#' evaluateMF(object=gaussian2, x=0)
#' derivateMF(object=gaussian2, x=0, i=1)
#' derivateMF(object=gaussian2, x=0, i="mu")
setMethod(f="derivateMF",signature = "GaussianMF", 
  definition=function(object,x,i){
  out<-switch(i,
    mu=1/(sqrt(2*pi)*object@parameters["sigma"]^3)*
      exp(-1/2*((x-object@parameters["mu"])/object@parameters["sigma"])^2)*
      (x-object@parameters["mu"]),
    sigma=(exp(-1/2*((x-object@parameters["mu"])/
      object@parameters["sigma"])^2)/sqrt(2*pi))*
      ((((x-object@parameters["mu"])^2)/object@parameters["sigma"]^4)-
      (1/object@parameters["sigma"]^2))
  )
  ##Check option
  if(is.null(out)){
    stop("Not correct index i for the partial derivate")
  }
  if(is.nan(out)){
    stop(paste("Partial derivate returned NaN for x=", x, " i=", i, 
      " parameters=", toString(object@parameters),sep=""))
  }
  return(out)
})
#'
#' @name derivateMF
#' @rdname derivateMF
#' @aliases derivateMF,NormalizedGaussianMF-method
#' @inheritParams derivateMF
#' @examples
#' #NormalizedGaussianMF example I
#' #A normalized Gaussian membership function with default parameters (mu=0, sigma=1)
#' #The derivate of the first parameter at x, should be 1
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' normalizedGaussian <- new(Class="NormalizedGaussianMF")
#' normalizedGaussian
#' evaluateMF(object=normalizedGaussian, x=0)
#' derivateMF(object=normalizedGaussian, x=0, i=1)
#' derivateMF(object=normalizedGaussian, x=0, i="mu")
#' #
#' #NormalizedGaussianMF example II
#' #A normalized Gaussian membership function with parameters (mu=0, sigma=1)
#' #The derivate of the first parameter at x, should be 1
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' normalizedGaussian2 <- new(Class="NormalizedGaussianMF",
#'  parameters=c(mu=0,sigma=1))
#' normalizedGaussian2
#' evaluateMF(object=normalizedGaussian2, x=0)
#' derivateMF(object=normalizedGaussian2, x=0, i=1)
#' derivateMF(object=normalizedGaussian2, x=0, i="mu")
setMethod(f="derivateMF",signature = "NormalizedGaussianMF", 
  definition = function(object,x,i){
  out <- switch(i,
    mu=1/(object@parameters["sigma"]^2)*
      exp(-1/2*((x-object@parameters["mu"])/object@parameters["sigma"])^2)*
      (x-object@parameters["mu"]),
    sigma=1/(object@parameters["sigma"]^3)*
      exp(-1/2*((x-object@parameters["mu"])/object@parameters["sigma"])^2)*
      (x-object@parameters["mu"])^2
  )
  ##Check option
  if(is.null(out)){
    stop("Not correct index i for the partial derivate")
  }
  if(is.nan(out)){
    stop(paste("Partial derivate returned NaN for x=", x, " i=", i, 
      " parameters=", toString(object@parameters),sep=""))
  }
  return(out)
})

