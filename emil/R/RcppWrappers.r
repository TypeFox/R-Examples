#' @useDynLib emil
#' @importFrom Rcpp evalCpp
NULL

#' Check if an object contains more than one unique value
#' 
#' @param x Vector or matrix.
#' @param na.rm Whether to ignore missing values.
#' @return Logical scalar that is \code{TRUE} if \code{x} contains more than one
#'   unique value and \code{FALSE} if not.
#' @return A logical scalar that is \code{TRUE} if \code{x} contains more than
#'   one unique value and \code{FALSE} otherwise. In case \code{x} contains
#'   missing values \code{NA} is returned if \code{na_rm = FALSE}. If there are
#'   no non-missing values \code{NA} is always returned.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
is_constant <- function(x, na.rm=FALSE){
    switch(base::mode(x),
        character = is_constant_character(x, na.rm),
        complex = is_constant_complex(x, na.rm),
        logical = is_constant_numeric(x, na.rm),
        numeric = is_constant_numeric(x, na.rm),
        stop(sprintf("Objects of mode `%s` are not supported.", base::mode(x))))
}

