#########################################
#Print method for DiagTest3Grp Object
#########################################
print.DiagTest3Grp <- function (x,...) 
{
  ####x: a DiagTest3Grp object
  cat("The DiagTest3Grp summary measure: ", x$type,"\n")
  cat("\n\n")
  cat(paste("Method used for ",x$type,":",x$method,"\n",sep=""))
  cat("\n\n")
  cat("Raw Data Summary:\n")
  print(x$dat.summary,...)
  cat("\n\n")
  cat(paste(x$type, "=", round(x$estimate, digits = 4), ",", (1-x$alpha)*100,"%", " CI=", round(x$CI[1], digits = 4),"~",round(x$CI[2], digits = 4),"\n",sep="") )
  cat("\n\n")
  cat(paste("Best cut-points: lower=", round(x$cut.point[1],digits=4),", upper=",round(x$cut.point[2],digits=4),sep=""),"\n")
  cat("\n\n")
  cat("The group correct classification probabilities are:\n")
  print(round(x$classify.prob,digits=4))      
  cat("\n\n")
  if(!is.na(x$sampleSize)) cat(paste("Sample Size to estimate ",x$type, " within specified margin of error=",x$sampleSize,"\n",sep=""))
}

