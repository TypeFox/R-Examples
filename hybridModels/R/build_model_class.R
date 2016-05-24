#' It builds an object of a pre-specified class.
#' 
#' @description  \code{buildModelClass} is generic function that calls a method
#'               to create a object base on model's name.
#' 
#' @param x is an empty object of a class requested.
#' 
#' @inheritParams hybridModel
#' 
#' @return An object of the class requested.
#' 
#' @export
#' @references .
buildModelClass <- function(x, var.names, init.cond, model.parms, prop.func = NULL,
                            state.var = NULL, state.change.matrix = NULL) UseMethod("buildModelClass")