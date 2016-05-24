#' Discriminant of the elliptic curve \eqn{y(x)}
#' 
#' Discriminant of the elliptic curve \eqn{y(x)}
#'
#' @method discr ecd
#'
#' @param object an object of ecd class
#' @param no.validate logical, if \code{TRUE}, don't validate presence of \code{beta}.
#'                    Default is \code{FALSE}.
#'
#' @return the discriminant
#'
#' @keywords stats
#'
#' @author Stephen H-T. Lihn
#'
#' @export discr
#'
#' @examples
#' d <- ecd(-1,1)
#' discr(d)
#'
### <======================================================================>
"discr.ecd" <- function(object, no.validate=FALSE)
{
    if(! no.validate) {
        if(object@beta != 0){
            stop("Parameter 'beta' must be zero for discriminant!\n")
        }
    }
    -16 * (4*object@gamma^3 + 27*object@alpha^2)
}
### <---------------------------------------------------------------------->
#' @rdname discr.ecd
setGeneric("discr", function(object, no.validate=FALSE) standardGeneric("discr"))
#' @rdname discr.ecd
setMethod("discr", signature("ecd"), discr.ecd)
### <---------------------------------------------------------------------->
