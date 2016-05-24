#' Print a parameterDef Object in the console
#' @name print.parameterDef
#' @aliases print.parameterDef
#' @title Print a parameterDef Object.
#' @method print parameterDef
#' @param x A parameterDef Object
#' @param \dots unused
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @seealso \code{\link{createParDef}}

print.parameterDef <-
function(x,...){
    cat("Selection Parameters:\n")
    print(x$selection)
    cat("\n")
    cat("Banker Parameters:\n")
    print(x$banker)
}
