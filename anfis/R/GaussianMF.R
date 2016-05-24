#' GaussianMF Membership Function S4 class
#'
#' Represent a concrete GaussianMF shaped Membership Function S4 class with 
#' parameters mu, sigma. Slots inherited of MembershipFunction class and related
#' functions: show, print, derivateMF, evaluateMF, [ and [<-.
#'
#' @slot parameters named numeric vector with parameters of Membership Function.
#' @slot nParameters integer with the number of parameters for validity check.
#' @slot name character The description of the membership function.
#' @slot expression expression object just to display purposes.
#'
#' @note derivateMF, evaluateMF are extended. Prototype is defined and validity 
#' is inherited.
#'
#' @include BellMF.R
#' @name GaussianMF-class
#' @rdname GaussianMF-class
#' @exportClass GaussianMF
#' @seealso \code{\link{BellMF-class}} and 
#'  \code{\link{NormalizedGaussianMF-class}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
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
#' #A Gaussian membership function with parameters (mu=0, sigma=1)
#' #The membership of x in the Gaussian, should be 1/sqrt(2*pi) = 0.3989423
#' #The derivate of the first parameter at x, should be 0
#' #The derivate on "mu" parameter at x, should be 0 
#' gaussian2 <- new(Class="GaussianMF",parameters=c(mu=0,sigma=1))
#' gaussian2
#' evaluateMF(object=gaussian2, x=0)
#' derivateMF(object=gaussian2, x=0, i=1)
#' derivateMF(object=gaussian2, x=0, i="mu")
GaussianMF<-setClass(Class="GaussianMF", contains="MembershipFunction",
  prototype=prototype(
    parameters=c(mu=0,sigma=1),
    nParameters=2L, 
    name="Gaussian Membership Function",
    expression=expression(1/sqrt(2*pi*sigma^2)*exp(-1/2*((x-mu)/sigma)^2)))
)
