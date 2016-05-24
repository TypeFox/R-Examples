#' Pretty Printing of Spectral's Output
#'
#' Pretty printing of spectral's output.
#' 
#' @param x an object of class spectral
#' @param digits digits to round to
#' @param ... additional parameters
#' @usage \method{print}{spectral}(x, digits = 3, ...)
#' @return Invisible string of the printed object.
#' @export print.spectral
#' @examples
#' 
#' # see ?spectral
#' 
#' 
print.spectral <- function(x, digits = 3, ...){
	
  ## print call
  cat("Call:\n")  
  print(x$call)
  cat("\n")
  
  ## print fit
  x$showFit()

  ## print stats and pvalues, etc.
  mat <- round(x$p.value, digits = digits)
  mat <- format(mat, digits = digits)  
  mat <- cbind("", rbind("", mat))
  mat[,1] <- c("", "P(samp)", "Pearson X^2", "Likelihood G^2", "Freeman-Tukey", "Cressie-Read", "Neyman X^2")
  mat[1,] <- colnames(mat)  
  mat[1,1] <- "P Values"
  mat <- apply(mat, 2, format, justify = "right")
  
  mat2 <- round(x$p.value.se, digits = digits)
  mat2 <- format(mat2, digits = digits)
  mat2 <- cbind("", rbind("", mat2))
  mat2[,1] <- row.names(mat2)
  mat2[1,] <- colnames(mat2)    
  mat2[1,1] <- "Std Errs"  
  mat2 <- apply(mat2, 2, format, justify = "right")  
  
  cat(
    apply(
      apply(cbind(mat, mat2), 2, format), 
      1, paste, collapse = "  "
    ), 
    sep = "\n"
  )

  invisible()
}
