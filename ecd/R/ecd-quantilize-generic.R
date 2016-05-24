#' Add the quantile data to the ecd object
#' 
#' Add the quantile data to the ecd object if it is not created yet.
#'
#' @method quantilize ecd
#'
#' @param object an object of ecd class
#' @param show.warning logical, if \code{TRUE}, display a warning message.
#'                    Default is \code{FALSE}.
#'
#' @return an object of ecd class with a newly generated ecdq object.
#'
#' @keywords distribution
#'
#' @author Stephen H-T. Lihn
#'
#' @export quantilize
#'
#' @examples
#' \dontrun{
#'    d <- ecd(-1,1)
#'    quantilize(d)
#' }
### <======================================================================>
"quantilize.ecd" <- function(object, show.warning=FALSE)
{
    if (! ecd.has_quantile(object)) {
        quantile(object) <- ecdq(object)
        if (show.warning) {
            warning("I am quantilizing ecd for you. This is slow. Consider ecd(with.quantile=T).")
        }
    }
    object
}
### <---------------------------------------------------------------------->
#' @rdname quantilize.ecd
setGeneric("quantilize", function(object, show.warning=FALSE) standardGeneric("quantilize"))
#' @rdname quantilize.ecd
setMethod("quantilize", signature("ecd"), quantilize.ecd)
### <---------------------------------------------------------------------->
