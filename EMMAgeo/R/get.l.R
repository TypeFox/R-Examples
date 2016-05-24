#' Generate a vector of weight transformation values from l_min to l_max.
#' 
#' This function generates a sequence of weight transformation values that
#' range from l_min (by default zero) to l_max (by default 95 % of the actual 
#' maximum possible value). It is a wrapper for the function 
#' \code{test.l.max()}.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param n Numeric scalar, length of the output vector (by default 10).
#' @param max Numeric scalar, fraction of the maximum value (by default 0.95).
#' @param min Numeric scalar, minimum value (by default zero).
#' @return Numeric vector of weight transformation values.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{test.l.max}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X, envir = environment())
#' 
#' ## infer l-vector
#' l <- get.l(X = X, n = 5, max = 0.8, min = 0.02)
#' 
#' @export get.l
get.l <- function(
  X, 
  n = 10, 
  max = 0.95,
  min = 0
){
  
  ## estimate maximum possible l-value
  l.max <- test.l.max(X = X, n = 10)
  
  ## calculate potentially smaller max value
  l.max <- l.max * max
  
  ## generate output sequence
  l <- seq(from = min, 
           to = l.max, 
           length.out = n)
  
  ## return output  
  return(l)
}