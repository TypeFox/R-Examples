setMethod("show",
    signature(object = "summary.cold"),
    function (object) 
    {
	cat("\nCall:\n")
	dput(object@call)
	coef <- object@coefficients
	nas <- is.na(coef[, 1])
	cnames <- names(coef[, 1][!nas])
	coef <- matrix(rep(coef[, 1][!nas], 4), ncol = 4)
	coef[, 1] <- 1:dim(coef)[[1]]
	coef[, 3] <- object@se[, 1][!nas]
	coef[, 4] <- round(coef[, 2]/coef[, 3], 3)
	dimnames(coef) <- list(cnames, c("Label", "Value", "Std. Error", "t value"))
	cat("\nNumber of profiles used in the fit: ", object@n.cases, "\n")
	cat("\nCoefficients:\t\n") 
	print(coef[ ,  ])
	cat("\nLog likelihood: ", round(object@log.likelihood, 4),"\n")
	cat("\nAIC: ", round(object@aic, 4),"\n")
	cat("\nMessage: ", object@message,"\n")
	})

