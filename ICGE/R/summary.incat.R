summary.incat <- function(object, ...){
# object must be a list with three components
   cat ("        INCA test    \n")
   k <- length(object$ProjectionsU)
   cat(" INCA statistic value =", object$StatisticW0,  "\n")
   cat("\n % of significative tests for alpha= " , object$alpha," : ", 
     format(object$Percentage_under_alpha, digits=3), "\n")
}