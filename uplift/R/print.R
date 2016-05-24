######################################################################
# Print method for class upliftRF
######################################################################

print.upliftRF <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Uplift random forest", "\n",sep="")
  cat("Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n", sep="")  
  cat("Split method: ", x$split_method, "\n", sep="")  
}

######################################################################
# Print method for class ccif
######################################################################

print.ccif <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Causal conditional inference forest", "\n",sep="")
  cat("Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n", sep="")  
  cat("Split method: ", x$split_method, "\n", sep="")  
}

