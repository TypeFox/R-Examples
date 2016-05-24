exact.cond=function(b,c,n,alpha){
		psidach=(b+c)/n
		tetadach=(b-c)/n
				exactci=function (x, n, conflev) 
				{
					alpha <- (1 - conflev)
					if (x == 0) {
						ll <- -1
						ul <- 1 - (alpha/2)^(1/n)
					}
					else if (x == n) {
						ll <- (alpha/2)^(1/n)
						ul <- 1
					}
					else {
						ll <- 1/(1 + (n - x + 1)/(x * qf(alpha/2, 2 * x, 2 * 
							(n - x + 1))))
						ul <- 1/(1 + (n - x)/((x + 1) * qf(1 - alpha/2, 2 * (x + 
							1), 2 * (n - x))))
					}
					c(ll, ul)
				}
		clopper=exactci(b,b+c,1-alpha)
		exact_low=(2*clopper[1]-1)*psidach
		exact_u=(2*clopper[2]-1)*psidach
		cint=c(exact_low,exact_u)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = tetadach)
		class(rval) <- "htest"
		return(rval)

}