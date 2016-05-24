#' @import methods
setGeneric("alk", function(object) standardGeneric("alk"))
setGeneric("N", function(object) standardGeneric("N"))
setGeneric("method", function(object) standardGeneric("method"))
setGeneric("parameters", function(object) standardGeneric("parameters"))
setGeneric("name", function(object) standardGeneric("name"))
setGeneric("description", function(object) standardGeneric("description"))
setMethod("alk", signature = "ALKr", function(object) object@alk)
setMethod("N", signature = "ALKr", function(object) object@N)
setMethod("method", signature = "ALKr", function(object) object@method)
setMethod("parameters", signature = "ALKr",
          function(object) object@parameters)
setMethod("name", signature = "ALKr", function(object) object@name)
setMethod("description", signature = "ALKr",
          function(object) object@description)

#' Extract elements of ALKr class
#' 
#' @name [
#' @aliases [,ALKr-method
#' @docType methods
#' @rdname extract-methods
setMethod("[", signature = "ALKr",
          function(x, i) {
            switch(i,
                   alk = x@alk,
                   N = x@N,
                   method = x@method,
                   parameters = x@parameters,
                   name = x@name,
                   description = x@description)})

# Setters
setGeneric("name<-", function(object, value) standardGeneric("name<-"))
setGeneric("description<-",
           function(object, value) standardGeneric("description<-"))
setReplaceMethod("name", signature = "ALKr",
                 function(object, value) {
                   object@name <- value
                   validObject(object)
                   return(object)})
setReplaceMethod("description", signature = "ALKr",
                 function(object, value) {
                   object@description <- value
                   validObject(object)
                   return(object)})

