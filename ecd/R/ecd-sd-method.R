#' Standard deviation, variance, mean, skewness, and kurtosis of ecd
#' 
#' Convenience wrappers around ecd's stats data
#'
#' @param object an object of ecd class
#'
#' @return numeric or mpfr
#'
#' @keywords stats
#'
#' @export ecd.sd
#' @export ecd.var
#' @export ecd.mean
#' @export ecd.skewness
#' @export ecd.kurt
#' @export ecd.kurtosis
#'
#' @examples
#' d <- ecd(-1,1)
#' ecd.sd(d)
#' ecd.var(d)
#' ecd.mean(d)
#' ecd.skewness(d)
#' ecd.kurt(d)
#'
### <======================================================================>
"ecd.sd" <- function(object)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    object@stats$stdev
}
### <---------------------------------------------------------------------->
#' @rdname ecd.sd
"ecd.var" <- function(object)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    object@stats$var
}
### <---------------------------------------------------------------------->
#' @rdname ecd.sd
"ecd.mean" <- function(object)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }
    
    object@stats$mean
}
### <---------------------------------------------------------------------->
#' @rdname ecd.sd
"ecd.skewness" <- function(object)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    object@stats$skewness
}
### <---------------------------------------------------------------------->
#' @rdname ecd.sd
"ecd.kurt" <- function(object)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }
    
    object@stats$kurtosis
}
### <---------------------------------------------------------------------->
#' @rdname ecd.sd
"ecd.kurtosis" <- function(object)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }
    
    object@stats$kurtosis
}
### <---------------------------------------------------------------------->

