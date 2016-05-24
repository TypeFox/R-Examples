print.Ord.logreg <-
function(x, ...)
{
 if(class(x)!="Ord.logreg")
   stop("x is not of class Ord.logreg")
 k<-length(x$model)
 cat("Number of classes: ", k+1)
 cat("\n")
 cat("Cross-validation: ", x$CV)
 cat("\n")
 cat("\n")
 cat("Ord.logreg model:")
 cat("\n")
 cat("\n")
 for (i in 1:k)
  {
  cat("Category", i,"Tree")
  cat("\n")
  print(x$model[[i]])
  cat("\n") 
  }
}
