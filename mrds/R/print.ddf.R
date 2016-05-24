#' Simple pretty printer for distance sampling analyses 
#'
#' Simply prints out summary of the model which was fitted. For more
#' detailed information see \code{\link{summary}}.
#'
#' @param x a \code{ddf} object 
#' @param ... not passed through, just for S3 compatibility.
#' @export
#' @aliases print.ddf
#'
#' @author David L. Miller
#' @export
print.ddf<-function(x, ...){

  cat("\nDistance sampling analysis object\n")
  print(summary(x))
  invisible()
}
