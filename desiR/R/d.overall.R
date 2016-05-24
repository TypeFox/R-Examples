#' @title Combine individual desirabilities
#'
#' @description Combines any number of desirability values into an overall
#' desirability.
#'
#' @details This function takes any number of individual desirabilities and
#' combines them with a weighted geometric mean to give an overall
#' desirability. The weights should be chosen to reflect the importance of the
#' variables. The values of the weights do not matter, only their relative
#' differences. Therefore weights of 4, 2, 1 are the same as 1, 0.5, 0.25. In
#' both cases the second weight is half of the first, and the third weight is a
#' quarter of the first.
#'
#' @param ... Any number of individual desirabilities.
#' @param weights Allows some desirabilities to count for more in the overall
#' calculation. Defaults to equal weighting.
#'
#' @return Numeric vector of desirability values.
#' 
#' @examples
#' set.seed(1)
#' x1 <- rnorm(1000, mean=100, sd =5) # generate data
#' x2 <- rnorm(1000, mean=100, sd =5) 
#' 
#' d1 <- d.high(x1, cut1=90, cut2=110, scale=1)
#' d2 <- d.low(x2, cut1=90, cut2=110, scale=1)
#'
#' D <- d.overall(d1, d2, weights=c(1, 0.5))
#' plot(rev(sort(D)), type="l")

d.overall <- function(..., weights = NULL){
  
  # get number of variables
  n <- length(list(...))

  # check whether any variables passed
  if(n == 0) stop("Some desirabilities must be included\n")

  # merge into matrix
  d.all <- cbind(...)

  if(min(d.all) < 0 | max(d.all) > 1) stop("Desirabilities must be between 0 and 1\n")

  # check lengths match
  if(!is.null(weights)) {
    if(ncol(d.all) != length(weights) ) {
      stop("Number of weights does not match number of desirabilities\n")
    }
  }
  
  # equal weights if none provided
  if(is.null(weights) == TRUE) {
    w <- rep(1,n)/n
  } else {
    w <- weights
  }
  
  # weighted geometric mean by row (omit missing)
  y <- apply(d.all, 1, function(x) exp(sum(w[!is.na(x)] * log(x[!is.na(x)])/
                                           sum(w[!is.na(x)])) )
             )
  
  return(y)
}
