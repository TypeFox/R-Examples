#' Whether the ecd object has quantile data or not
#' 
#' Whether the ecd object has quantile data or not. This is mostly for internal use.
#'
#' @param object an object of ecd class
#'
#' @return logical, whether the object has quantile data or not.
#'
#' @keywords distribution
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
### <======================================================================>
"ecd.has_quantile" <- function(object)
{
    return (! is.nan(object@quantile@N_seg))
}
### <---------------------------------------------------------------------->
