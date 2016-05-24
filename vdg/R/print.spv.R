#' Print Method for S3 \code{spv} classes 
#' 
#' Simple print methods for S3 classes \code{spv}, \code{spvlist}, \code{spvforlist} and \code{spvlistforlist}. See 
#' \code{\link{plot.spv}} for examples.
#' 
#' @aliases print.spv print.spvlist print.spvforlist print.spvlistforlist
#' @param x Object of class \code{spv} or \code{spvlist}
#' @param \dots Unimplemented
#' @author Pieter C. Schoonees
#' @export
#' @keywords print
print.spv <- function(x, ...){
  cat("\nObject of class 'spv'\n")
  cat("\nCall:\n",  paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Sample dimensions:\n",  nrow(x$sample), " columns and ", ncol(x$sample), " rows\n\n", sep = "")
  if(!is.null(as.list(x$call)$type)){ 
    if(as.list(x$call)$type %in% c("s", "S", "sphere")) stype <- "Spherical" 
    else stype <- "Cuboidal"
    cat("Design space type:\n", stype, "\n\n", sep = "")
  }
  cat("Summary of", ifelse(x$unscaled, "Unscaled Prediction Variance (UPV):\n", "Scaled Prediction Variance (SPV):\n"))
  print(summary(x$spv))
}