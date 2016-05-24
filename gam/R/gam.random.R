"gam.random" <-
function(x, y, w, df = sum(non.zero), sigma = 0)
{
	df.inv <- function(n, df, sigma = sum(n)/df - mean(n), iterations = 10
		)
	{
		if(df > length(n))
			return(0)
		current.df <- sum(n/(n + sigma))
		if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
			sigma
		else {
			sigma <- exp(log(sigma) + (current.df - df)/(sum((
				n * sigma)/(n + sigma)^2)))
			Recall(n, df, sigma, iterations - 1)
		}
	}
	nw <- tapply(w, x, sum)
	non.zero <- !is.na(nw)
	if(is.null(df))
		df <- sum(non.zero)
	if(sigma == 0)
		sigma <- df.inv(nw[non.zero], df)
	df <- sum(nw[non.zero]/(nw[non.zero] + sigma))
	fit <- tapply(w * y, x, sum)/(nw + sigma)
	var <- as.vector(w/(nw[x] + sigma))
	residuals <- as.vector(y - fit[x])
	list(x = seq(along = nw), y = fit, residuals = residuals, var = var,
		nl.df = df, sigma = sigma)
}
