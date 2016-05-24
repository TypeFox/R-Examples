#' Efficient binary search for character vectors
#' 
#' @param val  values
#' @param tab  table to find values in
#' @param L    lower bound
#' @param H    upper bound
#' 
#' @author Martin Morgan, Neal Fultz
#' @references \url{http://stackoverflow.com/questions/20133344/find-closest-value-in-a-vector-with-binary-search/} and
#'             \url{https://stat.ethz.ch/pipermail/r-help/2011-April/274182.html}
#' @export
#' @examples
#'  bsearch7(sample(letters, 5000, replace=TRUE), letters)
bsearch7 <-
     function(val, tab, L=1L, H=length(tab))
{
     n <- length(val)
     b <- matrix(c(L,H),n,2,byrow=TRUE) 
     i0 <- seq_along(val)

     repeat {
         M <- (b[,1] + b[,2]) %/% 2L
         i <- tab[M] > val
         b[i0 + i * n] <- M - i - i + 1L
         if(!any(b[, 2] >= b[, 1])) break;
     }
     b[,1] - 1L
}

