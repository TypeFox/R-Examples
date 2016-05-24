#' @method print miivs 
#' @export
print.miivs <- function(x,...){
  modeqns   <- x$df
  cat("Model Equation Information \n")
  cat("\n")
  print(modeqns, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
  cat("\n")
  cat("\n")
  if (x$miivs.out == TRUE){
    x <- x$eqns
    cat(paste("instruments <- '", sep=""), "\n")
      for(i in 1:length(x)){
        cat(paste(x[[i]]$DVobs, " ~ ", paste0(x[[i]]$IV, collapse = " + " ), sep=""), "\n")
      }
    cat(paste("'", sep=""), "\n")
  }
}