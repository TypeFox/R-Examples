print.incat <- function(x, ...){
# x must be a list with three components
   cat ("        INCA test    \n")
   k <- length(x$ProjectionsU)
   cat(" INCA statistic value =", x$StatisticW0,  "\n")
   cat("\n U projections values: \n")
   for (i in 1:k){
        cat("   U_", i," = ", x$ProjectionsU[i], " \n", sep="")
   }
   cat("\n % of significative tests for alpha= " , x$alpha," : ", 
     format(x$Percentage_under_alpha, digits=3), "\n")
}