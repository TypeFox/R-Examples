## ======================================================================
## utils
## ======================================================================

##' A convenience for declaring a default value for a field, in the
##' vein of \code{\link[methods]{prototype}} for S4 classes, except
##' the default value is quoted and evaluated lazily.
##'
##' @title Fields with prototypes
##' @param name The name of the field
##' @param class The class of the field
##' @param value Default value that when evaluated
##' initializes the field
##' @return A list suitable for use with \code{\link{setRefClass}}
##' @author Michael lawrence
##' @export
##' @examples
##' Brush.gen <- setRefClass("Brush",
##'                          fields = fieldWithPrototype("color", "character", "red"))
##' brush <- Brush.gen$new()
##' brush$color
##' brush$color <- "blue"
##' brush$color
fieldWithPrototype <- function(name, class, value) {
  .name <- paste(".", name, sep = "")
  .init <- paste(".init", name, sep = ".")
  value <- substitute(value)
  body <- substitute({
    if (missing(val)) {
      if (!length(.init)) {
        .name <<- value
        .init <<- TRUE
      }
      .name
    }
    else {
      if (!is(val, .class))
        stop("Cannot set an object of type '", class(val), "' on '", name,
             "', a field of type '", .class, "'")
      .name <<- val
      .init <<- TRUE
    }
  }, list(.name = as.name(.name), name = name, .class = class, value = value,
          .init = as.name(.init)))
  structure(list(as.function(c(alist(val=), body)), class, "logical"),
            names = c(name, .name, .init))
}

##' Declares a signal field that is lazily populated when the field is
##' first accessed. This avoids the need for the
##' constructor/initializer to explicitly create the signal.
##'
##' @title Declaring a signal field
##' @param expr The expression that names the signal and specifies its
##' signature. See the example.
##' @return A list of field definitions, suitable for passing to
##' \code{\link{setRefClass}}.
##' @author Michael Lawrence
##' @examples
##' setRefClass("Dataset", fields = c(elements = "list",
##'   declareSignal(elementsChanged(which))))
##' @export
declareSignal <- function(expr) {
  expr <- substitute(expr)
  name <- deparse(expr[[1]])
  expr[[1]] <- quote(Signal)
  do.call(fieldWithPrototype, list(name, "Signal", expr))
}


