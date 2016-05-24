#' @name scale0
#' @keywords scale
#' @aliases scaler
#' @author Sven E. Templer
#' @title Scale Numeric Values to Ranges
#' @description 
#' Scale numeric values to a range from 0 to 1 with the function
#' \code{scale0} or to a chosen range with \code{scaler}.
#' @param x Numeric vector to transform.
#' @param r Numeric vector of length 2 for range to scale values of 
#' \code{x} to.
#' @param b Numeric vector of length 2 to define the border of \code{x}
#' to use as scaling minimum and maximum.
#' @examples
#' #
#' 
#' scale0(0:10)
#' scale0(-1:3)
#' scale0(2:3)
#' 
#' scaler(0:10)
#' scaler(0:10, 1:2)
#' scaler(0:10, 1:2, c(0, 20))
#' 
#' #

#' @rdname scale0
#' @export scale0
scale0 <- function (x) {
  xmin <- min(x, na.rm=T)
  x <- x - xmin
  xmax <- max(x, na.rm=T)
  x <- x / xmax
  return(x)
}

#' @rdname scale0
#' @export scaler
scaler <- function (x, r = c(0, 1), b = range(x, na.rm = TRUE)) {
  rl <- r[1]
  ru <- r[2]
  bl <- b[1]
  bu <- b[2]
  if (bl > min(x, na.rm = TRUE))
    stop("Lower border in b is greater than minimum of x.")
  if (bu < max(x, na.rm = TRUE))
    stop("Upper border in b is less than maximum of x.")
  (ru - rl) * (x - bl) / (bu - bl) + rl
}

