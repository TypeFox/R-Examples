summary.autofitVariogram = function(object, ...)
# Provides a summary function for the autofitVariogram object
{
     cat("Experimental variogram:\n")
     print(object$exp_var, ...)
     cat("\nFitted variogram model:\n")
     print(object$var_model, ...)
	 cat("Sums of squares betw. var. model and sample var.")
	 print(object$sserr)
}
