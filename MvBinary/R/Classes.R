###################################################################################
##' Constructor of [\code{\linkS4class{MvBinaryResult}}] class
##'
##' This S4 class contains the results from the function \link{MvBinaryEstim}.
##'  
##' \describe{
##'   \item{alpha}{numeric. It indicates the marginal probability of that each variables are equal to 1.}
##'   \item{epsilon}{numeric. It indicates the dependency strength of each variables (between 0 and 1) with the other block variables.}
##'   \item{delta}{binary. It indicates the kind of dependency: two variables affiliated into the same block are positively correlated if their delta are equal and they are negatively correlated otherwise.}
##'   \item{blocks}{numeric. It indicates the block of each variable.}
##'   \item{nbparam}{numeric. It indicates the number of continuous parameters.}
##'   \item{loglike}{numeric. The model likelihood.}
##'   \item{bic}{numeric. The model BIC.}
##' }
##'
##' @examples
##'   getSlots("MvBinaryResult")
##'
##' @name MvBinaryResult-class
##' @rdname MvBinaryResult-class
##' @exportClass MvBinaryResult
##'

setClass(
  Class = "MvBinaryResult", 
  representation = representation(
    alpha="numeric",
    epsilon="numeric",
    delta="numeric",
    blocks="numeric",
    nbparam="numeric",
    loglike="numeric",
    bic="numeric"
  ), 
  prototype = prototype(
    alpha=numeric(),
    epsilon=numeric(),
    delta=numeric(),
    blocks=numeric(),
    nbparam=numeric(),
    loglike=numeric(),
    bic=numeric()
  )
)
