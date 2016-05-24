#' Bell Membership Function S4 class
#'
#' Represent a concrete Bell shaped Membership Function S4 class with parameters
#'  a, b, c. Slots inherited of MembershipFunction class and related functions:
#'  show, print, derivateMF, evaluateMF, [ and [<-.
#'
#' @slot parameters named numeric vector with parameters of Membership Function.
#' @slot nParameters integer with the number of parameters for validity check.
#' @slot name character The description of the membership function.
#' @slot expression expression object just to display purposes.
#'
#' @note derivateMF, evaluateMF are extended. Prototype is defined and validity 
#' is inherited.
#'
#' @include MembershipFunction-show.R
#' @name BellMF-class
#' @rdname BellMF-class
#' @exportClass BellMF
#' @seealso \code{\link{GaussianMF-class}} and 
#'  \code{\link{NormalizedGaussianMF-class}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
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
BellMF<-setClass(Class="BellMF", contains="MembershipFunction",
  prototype=prototype(
    parameters=c(a=1, b=1, c=0), 
    nParameters=3L, 
    name="Bell Membership Function",
    expression=expression(1/(1 + (((x - c)/a)^2)^(b^2))))
)