#' ALKr
#'
#' Every function used to calculate Age-Length Keys returns an \code{ALKr}
#' object.
#'
#' \describe{
#'   \item{alk}{A \eqn{i \times j} matrix with the probability of an individual
#'     of length \eqn{i} having age \eqn{j}, i.e. \eqn{P(j|i)}}
#'   \item{N}{A \eqn{i \times j} matrix with the estimated number of individuals
#'     of length \eqn{i} and age \eqn{j}}
#'   \item{method}{A string with the name of the algorithm used to calculate the
#'     ALK}
#'   \item{params}{A named list with any parameters needed by the algorithm}
#'   \item{name}{A string with a user-defined name for the ALKr object}
#'   \item{description}{A string with a user-defined description for the ALKr
#'     object}
#' }
#'  
#' @name ALKr-class
#' @rdname ALKr-class
#' @exportClass ALKr
#' 
setClass("ALKr",
  representation(
    alk = "matrix",
    N = "matrix",
    method = "character",
    parameters = "list",
    name = "character",
    description = "character"),
  validity = function(object) {
    if (!identical(dim(object@alk), dim(object@N)))
      return("alk and N matrices have different dimensions")
    return(TRUE)
  }
)

setMethod("initialize", "ALKr",
  function(.Object, alk, N, method, parameters, name = "", description = "") {
    .Object@alk <- alk
    .Object@N <- N
    .Object@method <- method
    .Object@parameters <- parameters
    .Object@name <- name
    .Object@description <- description
    validObject(.Object)
    .Object
  }
)

setMethod("show",
  signature(object = "ALKr"),
  function(object) {
    if (object@name != "") cat(paste(object@name,"\n"))
    if (object@description != "") cat(paste(object@description,"\n\n"))
    print(object@alk)
    cat(paste("\nMethod:", object@method,"\n"))
    if (length(object@parameters) > 0)
      print(matrix(object@parameters, dimnames = list(names(object@parameters), "Value")))
  }
)

