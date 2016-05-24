##' The functions
##' 
##' \code{all}, \code{any},
##' 
##' \code{sum}, \code{prod},
##' 
##' \code{min}, \code{max}, 
##' 
##' \code{range}, and
##'
##' \code{is.na}
##' 
##' for \code{hyperSpec} objects.
##' 
##' All these functions work on the spectra matrix.
##' @name Summary
##' @docType methods
##' @rdname summary
##' @aliases Summary,hyperSpec-method Summary all,hyperSpec-method
##'   any,hyperSpec-method sum,hyperSpec-method prod,hyperSpec-method
##'   min,hyperSpec-method max,hyperSpec-method range,hyperSpec-method
##' @param x hyperSpec object
##' @param ... further objects
##' @param na.rm logical indicating whether missing values should be removed
##' @return \code{sum}, \code{prod}, \code{min}, \code{max}, and \code{range} return  a numeric,
##' \code{all}, \code{any}, and \code{is.na} a logical.
##' @seealso \code{\link[base]{Summary}} for the base summary functions.
##' @export
##' @examples
##' 
##' 	range (flu) 
##' 
setMethod ("Summary", signature (x = "hyperSpec"),
           function (x, ..., na.rm = FALSE){
             validObject (x)

             if ((.Generic == "prod") || (.Generic == "sum"))
               warning (paste ("Do you really want to use", .Generic, "on a hyperSpec object?"))

             ## dispatch also on the objects in ...
             x <- sapply (list (x[[]], ...), .Generic, na.rm = na.rm)

             callGeneric (x, na.rm = na.rm)
           }
           )

##' @rdname summary
##' @aliases is.na,hyperSpec-method
##' @seealso \code{\link[base]{all.equal}} and \code{\link[base]{isTRUE}}
##' @export
##' @examples
##' 
##' is.na (flu [,, 405 ~ 410]);
setMethod ("is.na", signature (x = "hyperSpec"),
           function (x) {
             is.na (x@data$spc)
           })
           

