# Verteilungsfunktion für die eindimensionale Randdichte f(x_n) 
# einer Truncated Multivariate Student t Distribution,
# by integrating out (n-1) dimensions.
#
# @param xn Vektor der Länge l von Punkten, an dem die Verteilungsfunktion ausgewertet wird
# @param i Index  (1..n) dessen Randdichte berechnet werden soll
# @param mean (nx1) Mittelwertvektor
# @param sigma (nxn)-Kovarianzmatrix
# @param df degrees of freedom parameter
# @param lower,upper Trunkierungsvektor lower <= x <= upper
ptmvt.marginal <- function(xn, n=1, mean=rep(0, nrow(sigma)), sigma=diag(length(mean)), df = 1, 
		lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)))
{
	# check of standard tmvnorm arguments
	cargs <- checkTmvArgs(mean, sigma, lower, upper)
	mean  <- cargs$mean
	sigma <- cargs$sigma
	lower <- cargs$lower
	upper <- cargs$upper
	
	if (n < 1 || n > length(mean) || !is.numeric(n) || length(n) > 1 ||  !n %in% 1:length(mean))
	{
		stop("n must be a integer scalar in 1..length(mean)")
	}
	
	# Anzahl der Dimensionen                
	k = length(mean)
	
	Fx     = numeric(length(xn))
	upper2 = upper
	alpha  = pmvt(lower = lower, upper = upper, delta = mean, sigma = sigma, df = df)
	for (i in 1:length(xn))
	{
		upper2[n] = xn[i]
		Fx[i]     = pmvt(lower=lower, upper=upper2, delta=mean, sigma=sigma, df = df)
	}
	return (Fx/alpha)
}