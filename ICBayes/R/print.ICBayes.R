print.ICBayes <-
function(x, ...){
   cat("Call:\n")
   print(x$call)
   cat("\ncoef #Estimted regression coefficients\n")
   print(x$coef)
   cat("\ncoef_ssd #Sample standard deviation\n")
   print(x$coef_ssd)
   cat("\ncoef_ci #Credible interval of coefficients\n")
   print(x$coef_ci)
   cat("\nLPML #Log pseudomarginal likelihood\n")
   print(x$LPML)
}
