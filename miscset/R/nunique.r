#' @name nunique
#' @keywords unique
#' @author Sven E. Templer
#' @title Amount and Index of Unique Values
#' @description
#' Return the index or amount of unique values in a vector.
#' @param x Numeric vector to transform.
#' @param na Logical, \code{TRUE} to include/count \code{NA}.
#' @param ... Arguments forwarded to \link{unique}.
#' @examples
#' #
#' 
#' v <- c("a","b","a", NA)
#' nunique(v)
#' nunique(v, FALSE)
#' uniquei(v)
#' uniquei(v, FALSE)
#' 
#' #

#' @rdname nunique
#' @export
nunique <- function (x, na = TRUE, ...) {
    
  if (is.factor(x)) {
    if (na)
      nlevels(x) + any(is.na(x))
    else
      nlevels(x)
  } else {
    if (na)
      length(unique(x, ...))
    else
      length(unique(x[!is.na(x)], ...))
  }
  
}

#' @rdname nunique
#' @export
uniquei <- function (x, na = TRUE, ...) {
  
  xu <- if (na) unique(x, ...) else unique(x[!is.na(x)], ...)
  match(xu, x)
  
}
