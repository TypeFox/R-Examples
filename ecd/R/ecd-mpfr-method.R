#' Wrapper to convert numeric to mpfr
#'
#' Convert numeric to mpfr for ecd calculations.
#' \code{ecd.mp1} is the constant 1 wrapped in mpfr class.
#' \code{ecd.erf} is a wrapper on \code{Rmpfr::erf}.
#' \code{ecd.devel} is a developer tool to size down intensive mpfr tests for CRAN. Set \code{ecd_devel} in R options or OS env to change its value.
#'
#' @param x a numeric vector or list. If \code{x} is mpfr class,
#'          it will be passed through.
#' @param precBits an integer for mpfr precBits.
#'                 Default is from \code{getOption("ecd.precBits")}.
#'
#' @return The \code{mpfr} object
#'
#' @keywords utility
#'
#' @export ecd.mpfr
#' @export ecd.mp1
#' @export ecd.erf
#' @export ecd.devel
#'
#' @importFrom Rmpfr mpfr
#' @importFrom Rmpfr erf
#' @importClassesFrom Rmpfr mpfr
#'
#' @examples
#' x <- ecd.mpfr(1)
#' y <- ecd.mpfr(c(1,2,3))
#' z <- ecd.mp1
### <======================================================================>
"ecd.mpfr" <- function(x, precBits=getOption("ecd.precBits"))
{
	if (is.null(precBits)) {
	    precBits <- 120L
	} else if (is.na(precBits) | is.nan(precBits) | precBits < 2L) {
	    precBits <- 120L
	}
	
    c <- class(x)
    if (c=="numeric") return(Rmpfr::mpfr(x, precBits))
    if (c=="mpfr") return(x)
    
    # the list is usually caused by mclapply (the apply family)
    if (c=="list") return(new("mpfr", unlist(x)))
    
    
    stop(paste("Unknown class to convert to mpfr:", c))
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.mp1 <- ecd.mpfr(1)
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.erf <- function(x) Rmpfr::erf(x)
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
# export ecd_devel=1 to reduce the scope when not in "devel" mode
# CRAN has a 20-minute limit on tests
ecd.devel <- function() {
    dev <- getOption("ecd_devel")
    if (! is.null(dev)) { if (dev == TRUE) return(TRUE) }
    if (nchar(Sys.getenv("ecd_devel"))>0) return(TRUE)
    return(FALSE)
}
### <---------------------------------------------------------------------->


