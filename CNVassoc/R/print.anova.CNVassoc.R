print.anova.CNVassoc<-function(x,...){
  cat("\n\n--- Likelihood ratio test comparing 2 CNVassoc models:\n\n")
  cat("\nModel 1 call: ", deparse(x$mod1$call), "\n")
  cat("\nModel 2 call: ", deparse(x$mod2$call), "\n\n")
  cat("Chi=",x$chi.test["chi"],"(df=",x$chi.test["df"],")"," p-value=",x$chi.test["pvalue"],"\n\n")
  cat("  Note: the 2 models must be nested, and this function doesn't check this!\n")
}
  