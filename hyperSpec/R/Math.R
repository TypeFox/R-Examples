##' Math Functions for hyperSpec Objects
##' 
##' The functions \code{abs}, \code{sign},
##' 
##' \code{sqrt},
##' 
##' \code{floor}, \code{ceiling}, \code{trunc}, \code{round}, \code{signif},
##' 
##' \code{exp}, \code{log}, \code{expm1}, \code{log1p},
##' 
##' \code{cos}, \code{sin}, \code{tan}, \code{acos}, \code{asin}, \code{atan},
##' \code{cosh}, \code{sinh}, \code{tanh}, \code{acosh}, \code{asinh},
##' \code{atanh},
##' 
##' \code{lgamma}, \code{gamma}, \code{digamma}, \code{trigamma},
##' 
##' \code{cumsum}, \code{cumprod}, \code{cummax}, \code{cummin}
##' 
##' for \code{hyperSpec} objects.
##' 
##' @aliases  Math Math2 Math,hyperSpec-method Math2,hyperSpec-method abs,hyperSpec-method
##' sign,hyperSpec-method sqrt,hyperSpec-method floor,hyperSpec-method ceiling,hyperSpec-method
##' trunc,hyperSpec-method round,hyperSpec-method signif,hyperSpec-method exp,hyperSpec-method
##' log,hyperSpec-method expm1,hyperSpec-method log1p,hyperSpec-method cos,hyperSpec-method
##' sin,hyperSpec-method tan,hyperSpec-method acos,hyperSpec-method asin,hyperSpec-method
##' atan,hyperSpec-method cosh,hyperSpec-method sinh,hyperSpec-method tanh,hyperSpec-method
##' acosh,hyperSpec-method asinh,hyperSpec-method atanh,hyperSpec-method lgamma,hyperSpec-method
##' gamma,hyperSpec-method digamma,hyperSpec-method trigamma,hyperSpec-method cumsum,hyperSpec-method
##' cumprod,hyperSpec-method cummax,hyperSpec-method cummin,hyperSpec-method round,hyperSpec-method
##' signif,hyperSpec-method
##' @param x the \code{hyperSpec} object
##' @param digits integer stating the rounding precision
##' @return a \code{hyperSpec} object
##' @rdname math
##' @author C. Beleites
##' @seealso \code{\link[methods]{S4groupGeneric}} for group generic methods.
##' 
##' \code{\link[base]{Math}} for the base math functions.
##'
##' \code{\link[hyperSpec]{Arith}} for arithmetic operators,
##'   \code{\link[hyperSpec]{Comparison}} for comparison operators, and
##'   \code{\link[hyperSpec]{Summary}} for group generic functions working on
##'   \code{hyperSpec} objects.
##' @keywords methods math
##' @export
##' @examples
##' 
##' 	log (flu) 
##' 

setMethod ("Math2", signature (x = "hyperSpec"),
           function (x, digits){
             validObject (x)
             
             x [[]] <- callGeneric (x[[]], digits)
             
             x
           }
           )

##' @rdname math
##' @param ... ignored
##' @param base base of logarithm
##' @export
##' @aliases log log,hyperSpec-method
setMethod ("log", signature (x = "hyperSpec"),
					 function (x, base = exp (1), ...){
					 	validObject (x)
					 	
					 	x [[]] <-  log (x[[]], base = base)
					 	x
					 }
)

##' @rdname math
##' @export
setMethod ("Math", signature (x = "hyperSpec"),
					 function (x){
					 	validObject (x)
					 	
					 	if (grepl ("^cum", .Generic) || grepl ("gamma$", .Generic))
					 		warning (paste ("Do you really want to use", .Generic, "on a hyperSpec object?"))
					 	
					 	x [[]] <- callGeneric (x[[]])
					 	x
					 }
)


