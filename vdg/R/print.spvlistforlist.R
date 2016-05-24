#' @rdname print.spv
#' @method print spvlistforlist
#' @export
print.spvlistforlist <- function(x, ...){
  fornms <- names(x)
  desnms <- names(x[[1]])
  cat("Object of class 'spvlistforlist'\n")
  cat("\nDesign names: \n", paste0(desnms, collapse = ", "), sep = "")
  cat("\n\nFormulae: \n", paste0(fornms, collapse = ", "), sep = "")
  cat("\n\nCall:\n", paste(deparse(x[[1]][[1]]$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\nSample dimensions:\n",  nrow(x[[1]][[1]]$sample), " columns and ", ncol(x[[1]][[1]]$sample), " rows\n\n", sep = "")
  cat("Summary of", ifelse(x[[1]][[1]]$unscaled, "Unscaled Prediction Variance (UPV):\n", "Scaled Prediction Variance (SPV):\n"))
  spv_lst <- lapply(x, function(y) do.call(cbind, lapply(y, "[[", "spv")))
  for(i in seq_along(fornms)){
    cat("\n\tFormula:", fornms[i], "\n")
    cat("\tDesigns:\n")
    print(summary(spv_lst[[i]]))
  }
}