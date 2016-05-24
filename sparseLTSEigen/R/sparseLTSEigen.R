# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## function to call RcppEigen back end
#' @rdname sparseLTSEigen-package
#' @export
#' @import Rcpp 
#' @import RcppEigen
#' @import robustHD
#' @useDynLib sparseLTSEigen
.CallSparseLTSEigen <- function(..., PACKAGE) {
  .Call(..., PACKAGE="sparseLTSEigen")
}
