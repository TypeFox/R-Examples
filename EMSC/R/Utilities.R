#' Orthogonal vectors
#' 
#' Create orthogonal set of vectors that resemble the original input vectors.
#' 
#' @param X a \code{matrix} with vectors as rows (default) or columns (see below).
#' @param dim an integer specifying which dimension is the object dimension.
#' @param re a logical indicating if vectors should be norm-scaled before
#' orthogonalization and rescaled afterwards (default = TRUE).
#' 
#' @details The input vectors are orthgonalized using singular value
#' decomposition. To make the resulting vectors similar to the input
#' vectors (not just any base for the same space) they are re-oriented
#' towards the original vectors using Procrustes rotations.
#' 
#' To force the procedure to handle vectors of unequal magnitudes similarilly
#' they are by default rescaled to norm vectors before orthogonalization and
#' rescaled afterwards. This can be overridden using the \code{re} paramter.
#' 
#' @seealso \code{\link{EMSC}} \code{\link{EMSC_model}}  \code{\link{plot.EMSC}}
#' 
#' @importFrom pracma procrustes
#' @export
orthogonalVectors <- function(X, dim = 1, re = TRUE){
  # Handle dimensions
  if(dim == 1){
    X <- t(X)
  }
  
  # Check inputs
  pn <- dim(X)
  if(pn[2] > pn[1]){
    stop(paste('More objects (', pn[2], ') than variables (', pn[1], ').', sep = ""))
  }
  
  # Pre-re-sizing to reduce dominiance due to norm?
  if(re){
    Xr <- X * rep(1/sqrt(colSums(X^2)), each = pn[1])
  } else {
    Xr <- X
  }
  
  # Orthogonalization
  u <- svd(Xr, nu = pn[2], nv = 1)$u
  
  # Re-orientation
  basis <- procrustes(Xr,u)$P
  
  # Re-sizing
  basis <- basis%*%diag(sqrt(colSums(X^2))/sqrt(colSums(basis^2)))
  
  # Handle dimensions
  if(dim == 1){
    basis <- t(basis)
  }
  basis
}
