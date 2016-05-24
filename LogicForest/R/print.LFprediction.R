print.LFprediction <-
function(x, ...)
{
 if(class(x)!="LFprediction")
    stop("x not of class LFprediction")
 prdct<-x$LFprediction
 prop<-x$proportion_one
 if(length(x)==2)
  {
  cat("OOB Predicted values\n")
  cat("\n")
  print.default(prdct, quote=FALSE)
  cat("\n")
  cat("Proportion of OOB trees that predict 1")
  cat("\n")
  print.default(prop, quote=FALSE)
  }
 if(length(x)==3)
  {
  cat("Predicted values\n")
  cat("\n")
  print.default(prdct, quote=FALSE)
  cat("\n")
  cat("Proportion of trees that predict 1")
  cat("\n")
  print.default(prop, quote=FALSE)
  }
}
