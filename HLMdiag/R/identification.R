# Identifying unusual points
#
# This function will identify the unusual points found when
# looking at the plot comparing the shrinkage estimates and
# the random coefficients found through OLS.
#
# @param formula a formula that can be used with \code{lm()}
# @param identify the percentage of points to identify as unusual
# @author Adam Loy \email{loyad01@@gmail.com}
identify_resid <- function(eb, ols, identify){
#	model <- lm(formula = formula)
	yhat <- eb
#	attr(yhat, "names") <- rownames(eb)
	res <- ols - yhat
#	res <- resid(model)
	res <- res[order(abs(res), decreasing = TRUE)]
	unusual <- rep(FALSE, length(res))
	n <- round(length(res) * identify)
	unusual[1:n] <- TRUE
	return(data.frame(ids = names(res), residual = res, unusual = unusual))
}
