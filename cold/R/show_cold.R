setMethod("show",
    signature(object = "cold"),
    function (object) 
    {
	cat("\nCall:\n") 
	dput(object@call)
	cat("\nNumber of profiles used in the fit:", object@n.cases,"\n")
	cat("\nNumber of Coefficients:", length(object@coefficients),"\n")
	cat("\nLog likelihood:", round(object@log.likelihood, 4),"\n\n")
	cat("Message: ", object@message,"\n")
	})

