#' @rdname print.spv
#' @method print spvlist
#' @export
print.spvlist <- function(x, ...){
  nms <- names(x)
  if(is.null(nms)) nms <- seq_along(x)
  cat("\nObject of class 'spvlist'\n")
  cat("\nDesign names: \n", paste0(nms, collapse = ", "), sep = "")
  cat("\n\nCall:\n", paste(deparse(x[[1]]$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\nSample dimensions:\n",  nrow(x[[1]]$sample), " columns and ", ncol(x[[1]]$sample), " rows\n\n", sep = "")
  if(!is.null(as.list(x$call)$type)){
    if(as.list(x[[1]]$call)$type %in% c("s", "S", "sphere")) stype <- "Spherical" 
    else stype <- "Cuboidal"
    cat("Design space type:\n", stype, "\n\n", sep = "")
  }
  cat("Summary of", ifelse(x[[1]]$unscaled, "Unscaled Prediction Variance (UPV):\n", "Scaled Prediction Variance (SPV):\n"))
  spv_df <- do.call(cbind, lapply(x, "[[", "spv"))
  colnames(spv_df) <- nms
  cat("\n\tDesigns:\n")
  print(summary(spv_df))
}