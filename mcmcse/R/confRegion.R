library(ellipse)
confRegion <- function(mcse.obj, which = c(1,2), level = .95)
{
	mat <- mcse.obj$cov
	n <- mcse.obj$nsim
	p <- 2
	b <- mcse.obj$size
	a <- floor(n/b)

	m <- ifelse(mcse.obj$method == "bm", a-1, n - b)
	crit <- ifelse(mcse.obj$large, qchisq(level, df = p)/n,
              exp(log(p) + log(m) - log(n) - log(m-p+1) + log(qf(level, p, m-p+1))) )

	mu <- mcse.obj$est

	return(ellipse(mat, centre = mu[which], t = sqrt(crit), which = which))       

}