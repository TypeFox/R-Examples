#' @export
print.parallelSVM <- function(x, ...){
  
  # Calculate average number of support Vectors
  nSV <- c()
  for (i in 1:length(x)){
    nSV <- c(nSV,x[[i]]$tot.nSV)
  }
  tot.nSV <- round(mean(nSV))
  
  # Convert the type
  tempie <- x[[1]]$type
  tempie <- ifelse(tempie == 0, "C-classficiation",
                   ifelse(tempie == 1,"nu-classification",
                          ifelse(tempie == 2, "one-classification",
                                 ifelse(tempie == 3, "eps-regression", "nu-regression"))))
  
  # Convert the kernel
  kernie <- x[[1]]$kernel
  kernie <- ifelse(kernie == 0, "linear",
                   ifelse(kernie == 1,"polynomial",
                          ifelse(kernie == 2, "radial", "sigmoid")))
  
  # Print the call, Parameters and average number of Support Vectors
  cat("\nCall:\n")
  cat(attributes(x)$call)
  
  cat("\n\n\nParameters: \n")
  cat("   SVM-Type:  ")
  cat(tempie)
  cat("\n SVM-Kernel:  ")
  cat(kernie)
  cat("\n       cost:  ")
  cat(x[[1]]$cost)
  cat("\n      gamma:  ")
  cat(x[[1]]$gamma)
  
  cat("\n\nAverage Number of Support Vectors: ")
  cat(tot.nSV)
  cat("\n ")
}