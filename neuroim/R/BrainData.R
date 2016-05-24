#' @include AllClass.R
{}
#' @include AllGeneric.R
{}


 
#' @export
#' @rdname ndim-methods
setMethod(f="ndim", signature=signature(x = "BrainData"),
          def=function(x) length(dim(x@space)))

  
#' dim of \code{BrainData} object
#' @param x the object
#' @export
setMethod(f="dim", signature=signature(x = "BrainData"),
          def=function(x) dim(x@space))

  

#' @export
#' @rdname space-methods
setMethod(f="space", signature=signature(x = "BrainData"),
          def=function(x) x@space)

  

#' @export
#' @rdname spacing-methods
setMethod(f="spacing",signature= signature(x = "BrainData"),
          def=function(x) {
            sp <- space(x)
            spacing(sp)
          })

#' convert \code{BrainData} instance to matrix
#' @param x the object
#' @export
setMethod(f="as.matrix", signature=signature(x = "BrainData"), def=function(x) as(x, "matrix"))

#' convert \code{BrainData} instance to array
#' @param x the object
#' @export
setMethod(f="as.array", signature=signature(x = "BrainData"), def=function(x) as(x, "array"))

#' convert \code{BrainData} instance to vector
#' @param x the object
#' @export
setMethod(f="as.vector", signature=signature(x = "BrainData"), def=function(x) as(x, "vector"))
