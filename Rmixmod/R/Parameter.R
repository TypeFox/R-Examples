###################################################################################
##                                  Parameter.R                                  ##
###################################################################################


###################################################################################
##' Constructor of [\code{\linkS4class{Parameter}}] class
##'
##' This class defines parameters of a Mixture Model.
##'
##' \describe{
##'   \item{proportions}{a numeric vector containing proportions of the mixture model.}
##' }
##'
##' @examples
##'   getSlots("Parameter")
##'
##' @name Parameter-class
##' @rdname Parameter-class
##' @exportClass Parameter
##'
setClass(
    Class="Parameter",
    representation=representation(
        proportions = "numeric",
        "VIRTUAL"
    ),
    prototype=prototype(
        proportions = numeric(0)
    )
)
###################################################################################
