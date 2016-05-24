#' Print an FPCA object
#'
#' Print a simple description of an FPCA object
#'
#' @param x An FPCA object.
#' @param ... Not used.
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' res
#'
#' @export
print.FPCA <- function(x, ...){
  obj = x;
  cat("Functional Principal Components Object for", tolower(obj$optns$dataType), "data.\n\n")
  cat("The optimal number of components selected is:", length(obj$lambda),"and \nthe first k (<=3) eigenvalues are: ");
  if ( length(obj$lambda) < 4) { 
    cat( round(obj$lambda,3) ,"\n");
  } else {
    cat( round(obj$lambda[1:3],3) ,"\n")
  }
}


 
