z.test <- function (x, 
			alternative = c("two.sided", "less", "greater"), 
    		mu = 0, sigma=1, conf.level = 0.95) 
{
    DNAME <- deparse(substitute(x))      # record name of data coming in
    alternative <- match.arg(alternative)      # fancy argument matching

	n <- length(x)
	x.bar <- mean(x, na.rm=T)
	se <- sigma/sqrt(n)

	alpha <- 1 - conf.level
	p <- switch(alternative, 
		two.sided = c( alpha/ 2, 1-alpha/2 ),
		less = c( 0, 1-alpha ),
		greater = c( alpha, 1 )
		)
	z.star <- qnorm( p )

    Z <- ( x.bar - mu ) / se; names(Z) <- "z"
    SIGMA <- sigma; names(SIGMA) <- "sigma"
    MU <- mu; names(MU) <- "mean"
    ESTIMATE <- x.bar; names(ESTIMATE) <- "sample mean"
    CINT <- x.bar + z.star * se; 
		attr(CINT, "conf.level") <- conf.level            
    PVAL <- switch(alternative,
		two.sided = 2 * pnorm( - abs(Z) ),
		less = pnorm( Z ),
		greater = 1 - pnorm( Z ) 
		)
		
    structure( list( statistic = Z, parameter = SIGMA, p.value = PVAL, 
        conf.int = CINT, estimate = ESTIMATE, null.value = MU, 
        alternative = alternative, method = "Z test for a mean", 
        data.name = DNAME),        
        class = "htest")
}
