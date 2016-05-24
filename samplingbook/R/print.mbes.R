print.mbes <-
function(x,...){          
 cat("\nmbes object: Model Based Estimation of Population Mean\n")
 cat("Population size N = ",x$info$N, ", sample size n = ",x$info$n,"\n\n",sep="")
 cat("Values for auxiliary variable: \n")
 for(i in 1:x$info$p){
   cat("X.mean.",i," = ",round(x$info$aux[i],4),", x.mean.",i," = ",round(x$info$x.mean[i],4),"\n",sep="")
 }
 if(x$call$method=="simple"|x$call$method=="all"){
   cat("----------------------------------------------------------------\n")
   cat("Simple Estimate\n\n", sep="" )
   cat("Mean estimate: ", round(x$simple$mean,4),"\n")
   cat("Standard error: ", round(x$simple$se,4),"\n\n")
   cat(100 * x$call$level, "% confidence interval [", x$simple$ci[1], ",",x$simple$ci[2], "]\n\n",sep="")
 }
 if(x$call$method=="diff"|x$call$method=="all"){
   cat("----------------------------------------------------------------\n")
   cat("Difference Estimate\n\n", sep="" )
   cat("Mean estimate: ", round(x$diff$mean,4),"\n")
   cat("Standard error: ", round(x$diff$se,4),"\n\n")
   cat(100 * x$call$level, "% confidence interval [", x$diff$ci[1], ",",x$diff$ci[2], "]\n\n",sep="")
 }
 if(x$call$method=="ratio"|x$call$method=="all"){
   cat("----------------------------------------------------------------\n")
   cat("Ratio Estimate\n\n", sep="" )
   cat("Mean estimate: ", round(x$ratio$mean,4),"\n")
   cat("Standard error: ", round(x$ratio$se,4),"\n\n")
   cat(100 * x$call$level, "% confidence interval [", x$ratio$ci[1], ",",x$ratio$ci[2], "]\n\n",sep="")
 }
 if(x$call$method=="regr"|x$call$method=="all"){
   cat("----------------------------------------------------------------\n")
   cat("Linear Regression Estimate\n\n", sep="" )
   cat("Mean estimate: ", round(x$regr$mean,4),"\n")
   cat("Standard error: ", round(x$regr$se,4),"\n\n")
   cat(100 * x$call$level, "% confidence interval [", x$regr$ci[1], ",",x$regr$ci[2], "]\n\n",sep="")
   cat("----------------------------------------------------------------\n")
   cat("Linear Regression Model:")
   print(summary(x$regr$model))
 }
 invisible(x)
}
