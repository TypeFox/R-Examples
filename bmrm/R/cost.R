

#' Compute or check the structure of a cost matrix 
#' 
#' @param y a factor representing the labels of the instances
#' @param C either a cost matrix to check for consistency with labels in y, or a character string defining the standard matrix to compute. 
#'        If a character string the accepted values are "0/1" for a 0-1 cost matrix or "linear" for linear cost.
#' @return the cost matrix object
#' @export
#' @seealso bmrm, ordinalRegressionLoss
costMatrix <- function(y,C=c("0/1","linear")) {
  y <- as.factor(y)
  if (is.character(C)) {
    C <- match.arg(C)
    C <- switch(C,
           "0/1" = {
              C <- matrix(1,nlevels(y),nlevels(y),dimnames = list(levels(y),levels(y)))
              diag(C) <- 0
              return(C)
           },
           "linear" = abs(outer(levels(y),levels(y),'-'))
    )
  } else {
    C <- as.matrix(C)
    if (nrow(C)!=ncol(C)) stop("C must be a square matrix")
    if (nlevels(y)!=nrow(C)) stop("dimension of the square matrix C doesn't match with number of levels in y")
  }
  return(C)
}

