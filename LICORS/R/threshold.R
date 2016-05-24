#' @title Threshold a matrix/vector below and above
#'
#' @description 
#' \code{threshold} sets values of a vector/matrix below \code{min} 
#' to \code{min}; values above \code{max} are set to \code{max}.
#' 
#' \code{threshold} is mainly used  to project sparsified weight vectors 
#' (\code{\link{sparsify_weights}}) back onto the 
#' probability simplex (thus \code{min = 0} and then \code{\link{normalize}}).
#'
#' @param x a numeric matrix(like object)
#' @param min minimum value
#' @param max maximum value
#' @keywords manip array
#' @export
#' @seealso \code{\link{normalize}}
#' @examples
#' print(threshold(c(1,4,2,-1,10), min = 0))

threshold <- function(x, min = -Inf, max = Inf){
  
  if (min > -Inf){
    x[x < min] <- min
  } 
  if (max < Inf){
    x[x > max] <- max
  }
  invisible(x)
}