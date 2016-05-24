#' Simple pretty printer for distance sampling analyses 
#'
#' Simply prints out a brief description of the model which was fitted. For more
#' detailed information use \code{\link{summary}}.
#'
#' @param x a distance sampling analysis (result from calling \code{\link{ds}}).
#' @param ... not passed through, just for S3 compatibility.
#' @aliases print.dsmodel
#'
#' @author David L. Miller
#' @export
print.dsmodel<-function(x, ...){

  cat("\nDistance sampling analysis object\n")
  cat("\nDetection function:\n",model.description(x$ddf),"\n")

  if(!is.null(x$ddf$Nhat)){
    cat("\nEstimated abundance in covered region:",x$ddf$Nhat,"\n")
  }

  invisible()
}
