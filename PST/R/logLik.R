## Retunrs the log likelihood of a VLMC model

setMethod("logLik", "PSTf", function(object) {

	n <- nrow(object@data)

	message(" [>] model fitted to ", n, " sequence(s) - ", nobs(object), " symbols")

	pstsum <- summary(object)
	res <- object@logLik

	class(res) <- "logLik"
	attr(res, "df") <- pstsum@freepar

	return(res)
}
)
 
 
