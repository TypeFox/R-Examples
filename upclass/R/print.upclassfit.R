.packageName <- 'upclass'

print.upclassfit <- function (x, ...) # not sure 
{
  if (attr(x, "class") != "upclassfit") 
    stop("Incorrect class")
  #cat("'", class(x)[1], "' model object:\n", sep = "")
  #M<-x$Best$modelName
  cat("Model Name: ",x[["Best"]]$modelName,"\n")

    if (!is.null(x[["Best"]]$test$tab)) {

      cat("Total Misclassified: ",x[["Best"]]$test$misclass,"\n")
      cat("Misclassification Rate:  ",round(x[["Best"]]$test$rate,digits=3),"%\n\n")
    }
 
 # print(names(x))
#  cat("\n\nAvailable Components:\n")
#  cat("\n_________________________________________\n")
#  cat("\n_________________________________________\n")
}
