#' Summary Method for Class "lsbclust"
#' 
#' Summarize a \code{lsbclust} object.
#' @param object An object of class 'lsbclust'.
#' @param digits The number of digits in the printed output.
#' @param \dots Unimplemented.
#' @method summary lsbclust
#' @export
summary.lsbclust <- function(object, digits = 3, ...){
  
  ## Determine slots
  ind <- !sapply(object[c("overall", "rows", "columns", "interactions")], is.null)
  type <- c("overall", "rows", "columns", "interactions")[ind]
  
  ## Print basic information
  cat("Summary of 'lsbclust' Object\n\n")
  cat("Model terms: ", paste(type, collapse = ", "), "\n\n")
  cat("delta:", object$delta, "\n\n")
  
  ## Total loss
  lossmat <- rbind(loss = object$loss, df = object$df)
  lossmat <- cbind(lossmat, total = rowSums(lossmat))
  cat("Total loss:", lossmat[1, "total"], "on", lossmat[2, "total"], "degrees-of-freedom\n\n")
  
  ## Print loss and degrees-of-freedom
  cat("Breakdown of unscaled loss and degrees-of-freedom:\n")
  print(t(lossmat), digits = min(getOption("digits"), digits))
  
  ## Message for further use
  cat("\nUse summary(object$interactions) for information on the interaction fit.")
}