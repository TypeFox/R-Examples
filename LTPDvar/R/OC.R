OC <-
function(p, n, k, type = c("exact", "napprox","ewmaSK","ewma2"),lam=1)
{
	type = match.arg(type)
	Lnapprox <- function(p, n, k) pnorm((qnorm(1 - p) - k)/((1/n) +
    	((k^2)/(2 * n - 2)))^0.5)
	Lt <- function(p, n, k) 1 - pt(q = k * n^0.5, df = n - 1,
    	ncp = qnorm(1 - p) * n^0.5)
	OC_emwa2=function(p,n,k,lam) {
	c4=(2/(n-1))^0.5*gamma(n/2)/gamma((n-1)/2);
	pnorm((qnorm(1-p) - k*c4)*	
	(1/((1/n)*lam/(2-lam)+k^2*(1-c4^2))
	)^0.5)}	
	if (type == "exact")
    	return(Lt(p, n, k))
	if (type == "napprox")
    	return(Lnapprox(p, n, k))
	if (type == "ewmaSK")
  	return(pnorm((qnorm(1-p) - k)*(n/(lam/(2-lam)))^0.5))
	if (type == "ewma2")
	OC_emwa2(p,n,k,lam)
}
