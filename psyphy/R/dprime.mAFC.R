`dprime.mAFC` <-
function(Pc, m) {
	# m an integer > 2, number of choices
	# Pc - probability correct choice (unbiased observer)
	m <- as.integer(m)
	if (m < 2) stop("m must be an integer greater than 1")
	if (!is.integer(m)) stop("m must be an integer")
	if (Pc <= 0 || Pc >= 1) stop ("Pc must be in (0,1)")
	
	est.dp <- function(dp){	
		pr <- function(x) 
			dnorm(x-dp) * pnorm(x)^(m-1)#0-1)		
		Pc - integrate(pr, lower = -Inf, upper = Inf)$value 
		}		
	dp.res <- uniroot(est.dp, interval = c(-10,10))
	dp.res$root	
}

