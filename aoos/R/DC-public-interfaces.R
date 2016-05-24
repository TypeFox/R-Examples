setClass("public", contains = "VIRTUAL")
setClass("publicFunction", contains = c("public", "function"))
setClass("publicValue", contains = "publicFunction")
setClass("publicEnv", contains = c("public", "list"))
setClass("private", slots = c(.Data = "ANY"))

setGeneric("getPublicRepresentation", function(obj) obj)
setMethod("getPublicRepresentation", "publicEnv", function(obj) {
  obj@.Data[[1]]
})

#' @rdname defineClass
#' @export
setGeneric("private", function(x) new("private", x))

#' @rdname defineClass
#' @export
setMethod("private", "public", function(x) x)

#' Constructors for public members
#' 
#' These functions are used internally. You should not rely on them. Use \code{\link{public}} instead.
#' 
#' @param x a default value
#' @param validity an optional validity function for the set method. Returns TRUE or FALSE.
#' @param fun function definition
#' @param name name of member in refernece object
#' 
#' @rdname publicInterface
#' @export publicFunction
publicFunction <- function(fun) {
  new("publicFunction", .Data = fun)
}

#' @export
#' @rdname publicInterface
publicValue <- function(x = NULL, validity = function(x) TRUE) {
  force(x); force(validity)
  new("publicValue", .Data = function(value) {
    if(missing(value)) {
      return(x)
    } else {
      if(validity(value)) {
        x <<- value
        invisible(x)
      } else {
        stop("Invalid value!")
      }
    }
  })
}

#' @rdname publicInterface
setMethod("$", "publicEnv", function(x, name) {
  mc <- match.call()
  mc[[2]] <- quote(x@.Data[[1]])
  eval(mc)
})

#' @param x an object made public
#' @param validity function to check the validity of an object
#' 
#' @rdname defineClass
#' @export
setGeneric("public", function(x = NULL, validity = function(x) TRUE) {
  # A method for class 'environment' will coerce it's argument to an environment
  # The original class will be lost. I don't see why. Hence the if-statement.
  if(inherits(x, "environment")) {
    new("publicEnv", list(x))
  } else {
    publicValue(x, validity)
  }
})

#' @rdname defineClass
#' @export
setMethod("public", c(x = "function"), function(x, validity) {
  publicFunction(x)
})

#' @rdname defineClass
#' @export
setMethod("public", c(x = "private"), function(x, validity) {
  x
})

#' @rdname defineClass
#' @export
setMethod("public", c(x = "public"), function(x, validity) {
  x
})

