print.incaest <- function(x, ...){
# x must be a list with two components
   cat("INCA statistic value: \n")
   cat(x$Wvalue, " \n")
   cat("\n U projection values: \n")
   for (i in 1:length(x$Uvalue)){
       cat("U_",i, " = ", x$Uvalue[i], " \n", sep="")
   }
}