##' Construct a variable definition with (optional) additional metadata
##'
##' @param data an R vector of data to convert to the Crunch payload format.
##' See code{\link{toVariable}} for how R types are converted. If \code{data}
##' is not supplied, you may instead supply \code{values}, which will not be
##' converted in any way, nor will extra type information be supplied. Only
##' send \code{values} if you know what you're doing. You may also omit both
##' \code{data} and \code{values} to create an empty variable on the server
##' (all values will be system missing "No Data").
##' @param ... additional metadata attributes to send.
##' @return a \code{VariableDefinition} object, ready to POST to Crunch.
##' @export
##' @examples
##' VariableDefinition(rnorm(5), name="Some numbers",
##'     description="Generated pseudorandomly from the normal distribution")
##' VarDef(name="Integers", values=1:5, type="numeric",
##'     description="When creating variable definitions with 'values', you must
##'     specify 'type', and categorical variables will require 'categories'.")
##' @seealso \code{\link{toVariable}}
VariableDefinition <- function (data, ...) {
    out <- list(...)
    if (!missing(data)) {
        out <- updateList(toVariable(data), out)
    }
    class(out) <- "VariableDefinition"
    return(out)
}

##' @rdname VariableDefinition
##' @export
VarDef <- VariableDefinition

setOldClass("VariableDefinition")
