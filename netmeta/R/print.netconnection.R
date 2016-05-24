print.netconnection <- function(x,
                                digits = max(4, .Options$digits - 3), ...){
  
  meta:::chknumeric(digits, single = TRUE)
  
  
  matitle(x)
  ##
  cat(paste("Number of studies: k = ", x$k, "\n", sep=""))
  cat(paste("Number of treatments: n = ", x$n, "\n", sep=""))
  cat(paste("Number of pairwise comparisons: m = ", x$m, "\n", sep=""))
  ##
  cat("Number of networks: ", x$n.subnets, "\n\n", sep = "")
  
  cat("Distance matrix:\n")

  D <- round(x$D.matrix, digits = digits)
  D[is.infinite(D)] <- "."
  prmatrix(D, quote = FALSE, right = TRUE)
  
  invisible(NULL)
}
