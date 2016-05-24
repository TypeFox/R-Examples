"print.summary.ds.mixture"<-function (x,...){
# print summary data for a ds.mixture model objects
# 
# Argument
#  x	- object returned from summary.ds
#
# Value
#  Prints summary information
   cat("\nSummary for ds.mixture object \n")
     
   cat("Transect type          : ", x$ttype,"\n")
   cat("Mixture components     : ", x$mix.terms,"\n")
   cat("Number of observations : ", x$n,"\n")
   cat("Distance range         :  0 -",x$width,"\n")
   cat("AIC                    : ", x$aic, "\n")
   cat("\nDetection function parameters", "\n")
   print(x$coeff$pars)
   if(x$mix.terms>1){
      cat("\nMixture proportions (phi scale): ", "\n")
      print(x$coeff$mix.prop)
   }
   cat("\n")

   # plot the P_a and Nhat stuff
   parameters=data.frame(Estimate=c(x$average.p,x$Nhat))
   row.names(parameters)=c("Average p", "N in covered region")
   parameters$SE=c(x$average.p.se,x$Nhat.se)
   parameters$CV=parameters$SE/parameters$Estimate
   print(parameters)

   invisible()
}
