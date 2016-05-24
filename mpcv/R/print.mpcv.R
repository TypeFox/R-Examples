print.mpcv <- function(x, ...)
{
  cat(paste("Cpv: ", sprintf("%3d", x$CpV), "%\n", sep=""))
  cat(paste(" PS: ", sprintf("%3d", x$PS), "%", sep="")) 
  cat("  variable: ",x$PSvar , "\n") 
  cat(paste(" PD: ", sprintf("%3d", x$PD), "%", sep=""))
  cat("  variable: ", x$PDvar, "\n") 
  
}
