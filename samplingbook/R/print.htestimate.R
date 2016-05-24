print.htestimate <-
function(x,...){
  cat("\nhtestimate object: Estimator for samples with probabilities proportional to size\n")
  if(x$call$method=="yg"){cat("Method of Yates and Grundy:\n\n")} 
  if(x$call$method=="ht"){cat("Method of Horvitz-Thompson:\n\n")}
  if(x$call$method=="ha"){cat("Method of Hajek (approximate variance):\n\n") } 
  if(x$call$method=="hh"){cat("Method of Hansen-Hurwitz (approximate variance):\n\n") }
  cat("Mean estimator: ",x$mean,"\n",sep="")
  cat("Standard Error: ",x$se,"\n\n",sep="")
  invisible(x)
}
