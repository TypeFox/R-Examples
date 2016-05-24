.packageName <- 'upclass'

summary.upclassfit <- function (object, ...) 
{
  x<-object
  if (attr(x, "class") != "upclassfit") 
    stop("Incorrect class")
  cat("Model Name\n",x[["Best"]]$modelName,"\n")
  
   cat("Log Likelihood\n",x[["Best"]]$ll,"\n")
  cat("Dimension\n",x[["Best"]]$d,"\n")
  cat("Ntrain\n",x[["Best"]]$Ntrain,"\n")  
   cat("Ntest\n",x[["Best"]]$Ntest,"\n")
  cat("bic\n",x[["Best"]]$bic,"\n")
  
    if (!is.null(x[["Best"]]$test$tab)) {

      cat("Total Misclassified: ",x[["Best"]]$test$misclass,"\n")
      cat("Misclassification Rate:  ",round(x[["Best"]]$test$rate,digits=3),"%\n")
    }

}
