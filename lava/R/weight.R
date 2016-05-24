##' @export
`Weight` <- function(x,...) UseMethod("Weight")

##' @export
Weight.default <- function(x,...) eval(x$weight)
