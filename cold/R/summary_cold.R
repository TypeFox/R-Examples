setMethod("summary",
    signature(object = "cold"),
    function(object,cov=FALSE,cor=FALSE)
	{
	cat("\nCall:\n")
	print(object@call)
	cat("\nNumber of profiles in the dataset: ", object@ni.cases, "\n")
	cat("\nNumber of profiles used in the fit: ", object@n.cases, "\n")
	cat("\nLog likelihood: ", round(object@log.likelihood, 4),"\n")
	cat("\nAIC: ", round(object@aic, 4),"\n")
	coef <- object@coefficients
	nas <- is.na(coef[, 1])
	cnames <- names(coef[, 1][!nas])
	coef <- matrix(rep(coef[, 1][!nas], 5), ncol = 5)
	coef.aux<-matrix(rep(coef[, 1][!nas], 5), ncol = 5)
	coef[, 1] <- 1:dim(coef)[[1]]
	coef[, 3] <- object@se[, 1][!nas]
	coef[, 4] <- round(coef[, 2]/coef[, 3], 3)
	coef.aux[,1]<-pnorm(coef[, 2]/coef[, 3])
	coef.aux[,2]<-1-pnorm(coef[, 2]/coef[, 3])
	for (i in 1:dim(coef)[[1]])
	{coef.aux[i,3]<-min(coef.aux[i, 1],coef.aux[i, 2])}
	coef[, 5] <- round(2*coef.aux[,3],6)
	dimnames(coef) <- list(cnames, c("Label", "Value", "Std. Error", "t value", "p-value"))
	cat("\nCoefficients:\t\n") 
	if(all(is.na(match(cnames, "omega"))))
	print(coef[ ,  ])
	else {print(coef[-dim(coef)[[1]],  ])
	cat("\nRandom effect (omega):\t\n") 
	print(coef[dim(coef)[[1]],2:(dim(coef)[[2]]-2) ])}
		if (cov){
		cat("\nCovariance of Coefficients: \n")
		print(object@covariance, digits = 2)}
		if (cor){
		cat("\nCorrelation of Coefficients: \n")
		print(object@correlation, digits = 2)}
	cat("\nMessage: ", object@message,"\n")
	})

