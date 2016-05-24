#' @encoding UTF-8
#' @title Round Numbers Without Leading Zeros
#' @description Given a numeric vector, round numbers with no leading
#' zeros. Something nice for a plot or publicatio.
#' @param x A numeric vector.
#' @param digits An integer for the number of digits to round to.
#' @param add Logical, whether additional digits are to be added if no number appears in the pre-set digit level, default is \code{FALSE}.
#' @param max The Maximum number of digits to be shown, only affects if \code{add=TRUE}.
#' @return A vector of the same length of \code{x}, but stored as string.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @export
#' @examples
#'  x = seq(0, 1, by=.1)
#' rounded(x)
#'
`rounded` <- function(x, digits=2, add=FALSE, max=(digits+2)){
  y <- round(x, digits=digits)
  yk <- format(y, nsmall=digits)
  nzero <- sum(unlist(y)==0)
  if(add==TRUE){
    while(nzero>0){
      zeros <- y==0
      digits <- digits+1
      y[zeros] <- round(x, digits=digits)[zeros]
      yk[zeros] <- format(y[zeros], nsmall=digits)
      nzero <- sum(y==0)
      if(digits>(max-1))
        nzero <- 0
    }
  }##--end of add zeros
  z <- sub("^([-]?)0[.]","\\1.", gsub(" +", "", yk))
  return(noquote(z))
}##--end of rounded
NULL
