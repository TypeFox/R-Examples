#' NormalizedGaussianMF Membership Function S4 class
#'
#' Represent a concrete NormalizedGaussianMF shaped [0,1] Membership Function S4
#' class with parameters mu, sigma. Slots inherited of MembershipFunction class 
#' and related functions: show, print, derivateMF, evaluateMF, [ and [<-.
#'
#' @slot parameters named numeric vector with parameters of Membership Function.
#' @slot nParameters integer with the number of parameters for validity check.
#' @slot name character The description of the membership function.
#' @slot expression expression object just to display purposes.
#'
#' @note derivateMF, evaluateMF are extended. Prototype is defined and validity 
#' is inherited.
#'
#' @include GaussianMF.R
#' @name NormalizedGaussianMF-class
#' @rdname NormalizedGaussianMF-class
#' @exportClass NormalizedGaussianMF
#' @seealso \code{\link{BellMF-class}} and \code{\link{GaussianMF-class}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
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
NormalizedGaussianMF<-setClass(Class = "NormalizedGaussianMF",
  contains="MembershipFunction",
  prototype=prototype(
    parameters=c(mu=0,sigma=1),
    nParameters=2L, 
    name="Normalized Gaussian Membership Function",
    expression=expression(exp(-1/2*((x-mu)/sigma)^2)))
)
