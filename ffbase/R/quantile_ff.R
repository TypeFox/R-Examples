#' quantiles 
#'
#' The function quantile produces quantiles corresponding to the given 
#' probabilities. The smallest observation corresponds to a probability of 0 and the largest to a probability of 1.
#' Current implementation doesn't use the \code{type} parameter of \code{\link{quantile}}. For large \code{ff} vectors the
#' difference between the types is (very) small. If \code{x} has been \code{\link{ffordered}}, quantile is fast, otherwise it is $n log(n)$.
#' @method quantile ff 
#' @param x \code{ff} vector
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param na.rm logical; if true, any NA and NaN's are removed from x before the quantiles are computed.
#' @param names logical; if true, the result has a names attribute. Set to FALSE for speedup with many probs.
#' @param ... currently not used
#' @export
#' @export quantile.ff
#' @importFrom stats quantile
quantile.ff <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, ...){
  N <- length(x)
  
  nms <- if (names) paste(100*probs, "%", sep="") 
         NULL
  
  qnt <- 1L + as.integer(probs * (N-1))
  #print(qnt)
  
  idx <- ffordered(x)
  
  ql <- x[idx[qnt]]
  names(ql) <- nms
  ql
}

# x <- ff(1000000:1)
# #x <- addffIndex(x)
# 
# quantile(x)
