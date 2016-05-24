###################################################################################
##                                   Model.R                                     ##
###################################################################################

###################################################################################
##' @include global.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{Model}}] class
##'
##' This class defines the Mixmod models.
##'
##' \describe{
##'   \item{listModels}{character containing a list of models.}
##'   \item{free.proportions}{logical to include models with free proportions. Default is TRUE.}
##'   \item{equal.proportions}{logical to include models with equal proportions. Default is FALSE.}
##' }
##'
##' @examples
##'   getSlots("Model")
##'
##' @name Model-class
##' @rdname Model-class
##' @exportClass Model
##'
setClass(
    Class="Model",
    representation=representation(
        listModels = "character",
        free.proportions = "logical",
        equal.proportions = "logical",
        "VIRTUAL"
    ),
    prototype=prototype(
        listModels = character(0),
        free.proportions = logical(0),
        equal.proportions = logical(0)
    )
)
###################################################################################


###################################################################################
##' @rdname print-methods
##' @aliases print print,Model-method
##'
setMethod(
  f="print",
  signature=c("Model"),
  function(x,...){
    cat("****************************************\n")
    cat("*** MIXMOD Models:\n")
    cat("* list = ", x@listModels, "\n")
    if ( x@free.proportions  & x@equal.proportions )
      cat("* This list includes models with free and equal proportions.\n")
    else if ( x@free.proportions  & !x@equal.proportions )
      cat("* This list includes only models with free proportions.\n")
    else if ( !x@free.proportions  & x@equal.proportions )
      cat("* This list includes only models with equal proportions.\n")
    cat("****************************************\n")
  }
)
###################################################################################


###################################################################################
##' @rdname show-methods
##' @aliases show show,Model-method
##'
setMethod(
  f="show",
  signature=c("Model"),
  function(object){
    cat("****************************************\n")
    cat("*** MIXMOD Models:\n")
    cat("* list = ", object@listModels, "\n")
    if ( object@free.proportions  & object@equal.proportions )
    cat("* This list includes models with free and equal proportions.\n")
    else if ( object@free.proportions  & !object@equal.proportions )
    cat("* This list includes only models with free proportions.\n")
    else if ( !object@free.proportions  & object@equal.proportions )
    cat("* This list includes only models with equal proportions.\n")
    cat("****************************************\n")
  }
)
###################################################################################

