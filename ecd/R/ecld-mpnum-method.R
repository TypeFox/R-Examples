#' Wrappers for ecld to maintain consistent type between mpfr and numeric
#'
#' Primarily to make sure \code{x} is converted to mpfr vector if it is not,
#' when \code{use.mpfr} is set.
#'
#' @param object an object of ecd class
#' @param x a vector of numeric or mpfr.
#' @param test logical, test of \code{ifelse}.
#' @param yes return values for true elements of test
#' @param no  return values for false elements of test
#' @param FUN the function to be applied to each element of x
#' @param ... optional arguments to \code{FUN}
#'
#' @return a numeric or mpfr vector
#'
#' @keywords utility
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.mpnum
#' @export ecld.ifelse
#' @export ecld.sapply
#' @export ecld.mclapply
#'
### <======================================================================>
"ecld.mpnum" <- function(object, x)
{
    if (object@use.mpfr) ecd.mpfr(x) else x
}
### <---------------------------------------------------------------------->
#' @rdname ecld.mpnum
"ecld.ifelse" <- function(object, test, yes, no) {
    rs <- ifelse(test, ecld.mpnum(object, yes), ecld.mpnum(object, no))
    ecld.mpnum(object, rs)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.mpnum
"ecld.sapply" <- function(object, x, FUN, ...) {
    rs <- sapply(x, FUN, ...)
    ecld.mpnum(object, rs)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.mpnum
"ecld.mclapply" <- function(object, x, FUN, ...) {
    rs <- parallel::mclapply(x, FUN, ...)
    ecld.mpnum(object, rs)
}
### <---------------------------------------------------------------------->
