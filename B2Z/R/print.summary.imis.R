######################################
#It prints the posterior summaries   #
#for an object of class summary.imis #
######################################

print.summary.imis <- function(x, digits = max(options()$digits - 4, 3),...){
   cat("\nPosterior Summaries: \n")
   print(x$summary, digits = digits)
 
   cat("\nNote: GSD is the geometric standard deviation, i.e., GSD(x) = exp(sqrt(x))\n\n")

   cat("\nPosterior Covariance Matrix: \n")
   print(x$PostCovMat, digits = digits)

   cat("\nDeviance Information Criterion: ")

   cat("\npD: ")
   dput(round(x$pD,digits))

   cat("Dbar: ")
   dput(round(x$Dbar,digits))

   cat("DIC: ")
   dput(round(x$DIC,digits))


   cat("\n\nSampler used: IMIS ")
   cat("\n------------------\n")
   cat("\nEffective Sample Size (ESS):")
   dput(round(x$ESS,digits))

   cat("\nMaximum importance weight: ")
   dput(round(x$maxw,digits))

   cat("\nVariance of the rescaled importance weights: ")
   dput(round(x$Vhat,digits))

   cat("\nEntropy of importance weights relative to uniformity: ")
   dput(round(x$Uhat,digits))

   cat("\nExpected Fraction of Unique Points: ")
   dput(round(x$EFUP,digits))
 
   cat("\nExpected number of unique points after re-sampling: ")
   dput(round(x$Qhat,digits))
}

