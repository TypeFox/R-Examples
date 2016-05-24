#' Exporting to Re-export
#' 
#' @name e2r
#' @param x A Leontief decomposed Inter-Country Input Output table as created by decompr
#' @export
#' @import decompr
#' @examples 
#' # load the decompr package
#' library(decompr)
#' 
#' # load the example data set
#' data(leather)
#' 
#' # create a leontief decomposed data set
#' l <- decomp(inter,
#'             final,
#'             countries,
#'             industries,
#'             out)
#'  
#'  # apply the Exporting to Re-export
#'  e2r( l )


e2r <- function( x ) {
  
  # read attributes
  k      <- attr(x, "k")
  i      <- attr(x, "i")
  # rownam <- attr(x, "rownam")
  G <- length(k)
  N <- length(i)
  
  # transform back to 2dim x 2dim matrix
  x <- matrix(x[,5], nrow=G*N, byrow=TRUE)
  
  # remove exports to self
  f <- rowSums( x - diagonals::fatdiag(diagonals::fatdiag(x, steps=G), steps=G ) )
  
  # divide by own exports
  for (j in 1:N) {
    s <- seq( ((j-1)*N + 1), j*N )
    f[s] <- f[s] / sum(colSums(x[,s]))
  }
  
  f <- as.data.frame(f)
  f <- cbind(rep(k, each = N),
             rep(i, times = G),
             f)
  rownames(f) <- NULL  
  names(f) <- c("country", "sector", "e2r")
  
  return(f)
  
}

