#' Print a description of a density surface model variance object
#' 
#' This method only provides a short summary, use the 
#' \code{\link{summary.dsm.var}} method for information.
#' 
#' @param x a \code{dsm} variance object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @export
#' @author David L. Miller
#' @seealso \code{\link{summary.dsm.var}}
#' @keywords utility
print.dsm.var<-function(x,...){

  print(summary(x))

}
