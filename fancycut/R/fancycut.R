#' Like \code{cut}, turn a vector of numbers into a factor
#'
#' @param x           a numeric vector
#' @param intervals   a character vector of intervals
#' @param buckets     a character vector of levels for the new factor
#'   these have a 1-1 correspondence with \code{intervals}
#' @param na.bucket   what level should NA values be given?
#' @param unmatched.bucket
#'   what level should numbers not covered by an interval be given?
#' @param out.as.factor
#'   default is TRUE
#'   Should the resulting vector be a factor?
#'   If FALSE will return a character vector.
#'
#' @examples
#' fancycut(-10:10, c('(0,2]','(2,5)','[5,10]'), c('Small','Medium','Large'))
#' fancycut(-10:10, c('[0,0]','(0,2]','(2,5)','[5,10]'), c('Zero','Small','Medium','Large'))
#' @export
fancycut <- function(x, intervals, buckets = intervals,
                     na.bucket = NA, unmatched.bucket = NA,
                     out.as.factor = TRUE) {

  # Make sure that intervals and buckets are the same length
  l <- length(intervals)
  if(l != length(buckets)) {
    stop('FancyCut requires a 1-1 map from intervals to buckets')
  }

  # Make suer that x is numeric
  if (!is.numeric(x))
    stop("'x' must be numeric")


  out <- rep(NA, length(x))
  for(index in 1:l) {
    i <- intervals[index]
    b <- buckets[index]
    n <- nchar(i[1])
    left <- substr(i, 1, 1)
    right <- substr(i, n, n)
    bounds <- strsplit(substr(i, 2, n - 1), ",")
    upper <- as.numeric(bounds[[1]][2])
    lower <- as.numeric(bounds[[1]][1])

    mask <- rep(FALSE, length(x))
    if(left == '[' & right == ']') {mask <- x >= lower & x <= upper}
    if(left == '[' & right == ')') {mask <- x >= lower & x <  upper}
    if(left == '(' & right == ']') {mask <- x >  lower & x <= upper}
    if(left == '(' & right == ')') {mask <- x >  lower & x <  upper}

    out[mask] <- b
  }

  out[is.na(out)]  <- unmatched.bucket
  out[is.na(x)]    <- na.bucket

  levels <- unique(c(buckets, na.bucket, unmatched.bucket))

  if(out.as.factor) {
    return(factor(
      out,
      levels = levels,
      exclude = NULL
    ))
  } else {
    return(out)
  }
}

