
#' Return or set gradient attribute
#' 
#' gradient attribute isused by loss/risk function to return the gradient of
#' the function at a given point together with the function value
#' 
#' @name gradient
#' @rdname gradient
#' @aliases gradient<-
#' @title Return or set gradient attribute
#' @param x any R object
#' @param value new gradient value to set
#' @param ... additional paramters 
#' @return attr(x,"gradient")
#' @export
gradient <- function(x,...) UseMethod("gradient")

#' @rdname gradient
#' @export
gradient.default <- function(x,...) attr(x, "gradient")

#' @rdname gradient
#' @export
"gradient<-" <- function(x,...,value) UseMethod("gradient<-")

#' @rdname gradient
#' @export
"gradient<-.default" <- function(x,...,value) {attr(x, "gradient") <- value;x}

