summary.autoKrige = function(object, ...)
# Provides a summary function for the autoKrige object
{
     cat("krige_output:\n")
     print(summary(object$krige_output))
     cat("\nexp_var:\n")
     print(object$exp_var)
     cat("\nvar_model:\n")
     print(object$var_model)
	 cat("Sums of squares betw. var. model and sample var.")
	 print(object$sserr)
}
