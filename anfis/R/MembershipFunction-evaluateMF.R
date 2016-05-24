#' \code{evaluateMF} evaluate membership
#'
#' Evaluate de membership of x to the object MembershipFunction heirs.
#' 
#' @param object MembershipFunction class heirs
#' @param x numeric of the MembershipFunction to be evaluated
#'
#' @return 0 <= numeric <=1 with the obtained membership value
#'
#' @include NormalizedGaussianMF.R
#' @exportMethod evaluateMF
#' @docType methods
#' @name evaluateMF
#' @rdname evaluateMF
#' @aliases evaluateMF-methods
#' @seealso \code{\link{MembershipFunction-class}} and \code{\link{derivateMF}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setGeneric(name="evaluateMF", def=function(object, x){
  standardGeneric("evaluateMF")
})
#'
#' @name evaluateMF
#' @rdname evaluateMF
#' @aliases evaluateMF,MembershipFunction-method
#' @inheritParams evaluateMF
setMethod(f="evaluateMF", signature="MembershipFunction", 
  definition=function(object,x){
  stop("Only heir of MembershipFunction class can evaluateMF!!")
})
#'
#' @name evaluateMF
#' @rdname evaluateMF
#' @aliases evaluateMF,BellMF-method
#' @inheritParams evaluateMF
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
setMethod(f="evaluateMF",signature = "BellMF", 
  definition = function(object,x){
  return(1/(1 + (((x - object@parameters["c"])/
    object@parameters["a"])^2)^(object@parameters["b"]^2)))
})
#'
#' @name evaluateMF
#' @rdname evaluateMF
#' @aliases evaluateMF,GaussianMF-method
#' @inheritParams evaluateMF
#' @examples
#' #GaussianMF example I
#' #A Gaussian membership function with default prototype (mu=0, sigma=1)
#' #The membership of x in the Gaussian, should be 1/sqrt(2*pi) = 0.3989423
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' gaussian <- new(Class="GaussianMF")
#' gaussian
#' evaluateMF(object=gaussian, x=0)
#' derivateMF(object=gaussian, x=0, i=1)
#' derivateMF(object=gaussian, x=0, i="mu")
#' #
#' #GaussianMF example II
#' #A Gaussian membership function with paramateres (mu=0, sigma=1)
#' #The membership of x in the gaussian, should be 1/sqrt(2*pi) = 0.3989423
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' gaussian2 <- new(Class="GaussianMF",parameters=c(mu=0,sigma=1))
#' gaussian2
#' evaluateMF(object=gaussian2, x=0)
#' derivateMF(object=gaussian2, x=0, i=1)
#' derivateMF(object=gaussian2, x=0, i="mu")
setMethod(f="evaluateMF",signature = "GaussianMF", 
  definition = function(object,x){
  return(1/sqrt(2*pi*object@parameters["sigma"]^2)*
    exp(-1/2*((x-object@parameters["mu"])/object@parameters["sigma"])^2))
})
#'
#' @name evaluateMF
#' @rdname evaluateMF
#' @aliases evaluateMF,NormalizedGaussianMF-method
#' @inheritParams evaluateMF
#' @examples
#' #NormalizedGaussianMF example I
#' #A normalized Gaussian membership function with default paramateres (mu=0, sigma=1)
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
#' #A normalized Gaussian membership function with paramateres (mu=0, sigma=1)
#' #' #The derivate of the first parameter at x, should be 1
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0
#' normalizedGaussian2 <- new(Class="NormalizedGaussianMF",
#'  parameters=c(mu=0,sigma=1))
#' normalizedGaussian2
#' evaluateMF(object=normalizedGaussian2, x=0)
#' derivateMF(object=normalizedGaussian2, x=0, i=1)
#' derivateMF(object=normalizedGaussian2, x=0, i="mu")
setMethod(f="evaluateMF",signature = "NormalizedGaussianMF", 
  definition = function(object,x){
  return(exp(-1/2*((x-object@parameters["mu"])/object@parameters["sigma"])^2))
})
