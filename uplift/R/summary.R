######################################################################
# uplift method for summary 
######################################################################

summary.upliftRF <- function(object, ...) {
  
  res <- list (call         = object$call,
               importance   = varImportance(object, ...), 
               ntree        = object$ntree,
               mtry         = floor(object$mtry), 
               split_method = object$split_method)             
  
  class(res) <- "summary.upliftRF"
  return(res)
}


print.summary.upliftRF <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n", sep="")  
  cat("Split method: ", x$split_method, "\n", sep="")
  cat("\n")
  cat("Variable importance:\n")
  print(x$importance) 
}

######################################################################
# ccif method for summary 
######################################################################

summary.ccif <- function(object, ...) {
  
  res <- list (call         = object$call,
               importance   = varImportance(object, ...), 
               ntree        = object$ntree,
               mtry         = floor(object$mtry), 
               split_method = object$split_method)             
  
  class(res) <- "summary.ccif"
  return(res)
}


print.summary.ccif <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n", sep="")  
  cat("Split method: ", x$split_method, "\n", sep="")
  cat("\n")
  cat("Variable importance:\n")
  print(x$importance) 
}

