`print.mws` <-
function(x,...)
#print method for objects of class "mws"
#returns group membership for each session/subject
{
  cat("\n")
  if (x$method=="separate") {
    cat("Fitted model: Weibull mixture model \n")
  } else if (x$method=="propgroup") {
    cat("Fitted model: Weibull proportional hazards model (proportional group hazards) \n")
  } else if (x$method=="propsite") {
    cat("Fitted model: Weibull proportional hazards model (proportional variable hazards) \n")
  } else if (x$method=="propallmain") {
    cat("Fitted model: Weibull proportional hazards model (main effects) \n")
  } else if (x$method=="propallint") {
    cat("Fitted model: Weibull proportional hazards model (full model) \n")
  }
  cat("Number of model parameters:",x$npar,"\n")
  cat("\n")
  cat("Number of clustered objects:",length(x$group),"\n")
  cat("Number of clusters:",x$K,"\n")
  cat("\n")
  cat("Cluster means: \n")
  #colnames(x$clmean) <- paste("V",1:dim(x$clmean)[2],sep="")
  #rownames(x$clmean) <- paste("Cluster",1:x$K,sep="")
  print(x$clmean)
  cat("\n")
}

