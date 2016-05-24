`print.RRT`=function(x, ...){
  cat("\nCall:\n")
  print(x$Call)
  
  cat("\n",x$Model," model", sep = "")
  cat("\n",x$Name," model for the ",x$Type," estimator\nParameters: ",
       paste(names(x$Param),gsub(" ","",formatC(x$Param,digits=2)),sep="=",collapse="; "),"\n",sep="")

  if((x$Model=="Qualitative")&(x$Type=="total")){
    cat("\nEstimation:",round(x$Estimation))
  }else{
  cat("\nEstimation:",x$Estimation)
  }
  cat("\nVariance:",x$Variance)
  cat("\nConfidence interval (",format(100*(x$ConfidenceLevel),digits=2),"%)\n",sep="")
  cat("    Lower bound:",x$ConfidenceInterval[1],"\n")
  cat("    Upper bound:",x$ConfidenceInterval[2],"\n\n")
}