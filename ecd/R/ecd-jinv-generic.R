#' J-invariant of the elliptic curve \eqn{y(x)}
#' 
#' J-invariant of the elliptic curve \eqn{y(x)}
#'
#' @method jinv ecd
#'
#' @param object an object of ecd class
#' @param no.validate logical, if \code{TRUE}, don't validate presence of \code{beta}.
#'                    Default is \code{FALSE}.
#'
#' @return the j-invariant
#'
#' @keywords stats
#'
#' @author Stephen H-T. Lihn
#'
#' @export jinv
#'
#' @examples
#' d <- ecd(1,1)
#' j <- jinv(d)

### <======================================================================>
"jinv.ecd" <- function(object, no.validate=FALSE)
{
    if(! no.validate) {
        if(object@beta != 0){
            stop("Parameter 'beta' must be zero for j-inv!\n")
        }
    }
    
    D <- discr.ecd(object, no.validate=no.validate)
    ifelse(D==0, NaN, -1728 * (4*object@gamma)^3/D)
}
### <---------------------------------------------------------------------->
#' @rdname jinv.ecd
setGeneric("jinv", function(object, no.validate=FALSE) standardGeneric("jinv"))
#' @rdname jinv.ecd
setMethod("jinv", signature("ecd"), jinv.ecd)
### <---------------------------------------------------------------------->
