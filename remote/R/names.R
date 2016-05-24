if (!isGeneric('names')) {
  setGeneric('names', function(x, ...)
    standardGeneric('names')) 
}

#' Names of Eot* objects
#' 
#' @param x a EotMode or EotStack
#' 
#' @description
#' Get or set names of Eot* objects
#' 
#' @return
#' if \code{x} is a EotStack, the names of all mdoes, 
#' if \code{x} is a EotMode, the name the respective mode
#' 
#' @examples
#' data(vdendool)
#' 
#' nh_modes <- eot(vdendool, n = 2)
#' 
#' ## mode names
#' names(nh_modes)
#' names(nh_modes) <- c("vdendool1", "vdendool2")
#' 
#' names(nh_modes)
#' names(nh_modes[[2]])
#' 
#' @export
#' @name names
#' @rdname names
#' @aliases names,EotStack-method


setMethod('names', signature(x = 'EotStack'), 
          function(x) {
            x@names
          }
)

#' @export
#' @name names<-
#' @rdname names
#' @aliases names<-,EotStack-method

setMethod('names<-', signature(x = 'EotStack'), 
          function(x, value) {
            nm <- nmodes(x)
            if (is.null(value)) {
              value <- rep('', nm)
            } else if (length(value) != nm) {
              stop('incorrect number of mode names')
            }
            x@names <- value
            for (i in seq(x@modes)) {
              x@modes[[i]]@name <- x@names[i]
            }
            return(x)
          }
)

#' @export
#' @name names
#' @rdname names
#' @aliases names,EotMode-method

setMethod('names', signature(x = 'EotMode'), 
          function(x) { 
            x@name
          }
)

#' @export
#' @name names<-
#' @rdname names
#' @aliases names<-,EotMode-method
#' @param value name to be assigned

setMethod('names<-', signature(x = 'EotMode'), 
          function(x, value) {
            x@name <- value
            return(x)
          }
)
