if (!isGeneric('nmodes')) {
  setGeneric('nmodes', function(x, ...)
    standardGeneric('nmodes')) 
}

#' Number of modes of an EotStack
#' 
#' @param x an EotStack
#' 
#' @details
#' retrieves the number of modes of an EotStack
#' 
#' @return
#' integer
#' 
#' @examples
#' data(vdendool)
#' 
#' nh_modes <- eot(vdendool, n = 2)
#' 
#' nmodes(nh_modes)
#' 
#' @export
#' @name nmodes
#' @rdname nmodes
#' @aliases nmodes,EotStack-method

setMethod('nmodes', signature(x = 'EotStack'), 
          function(x) { 
            length(x@modes)
          }
)
