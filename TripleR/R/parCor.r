# partial correlation between x and y, correct for group membership z.
# z has to be a factor! (not a continuos variable)
# p values are corrected for loss of df
#' @export
parCor <- function(x,y,z) {
	if (sd(c(length(x), length(y), length(z))) != 0) stop("x, y, and z need to have the same length!")
	if (noVar(x) | noVar(y)) {
		warning("parCor: One or both variables have zero variance; skipping partial correlations!")
		return(NULL)
	}
	
	df <- na.omit(data.frame(x, y, z))
	
	z <- factor(as.character(z))	# reduce the factor to actually occuring levels
	x.z <- resid(lm(x~z, df, na.action=na.exclude))
	y.z <- resid(lm(y~z, df, na.action=na.exclude))
	
	par.cor <- cor(x.z, y.z, use="p")
	
	#Calculate the t value by hand:
	k <- length(levels(z)) - 1	# k = number of control parameters = number of groups - 1
	n <- nrow(df)				# n = number of participants
	dfs <- n-2-k
	
	t.value <- par.cor * sqrt(dfs/(1-(par.cor^2)))
	p.value <- (1-pt(abs(t.value), df=dfs, lower.tail = TRUE))*2   # always test two-sided

	return(list(par.cor=par.cor, t.value=t.value, df=dfs, p=p.value))
}