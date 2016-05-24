unif.test.quantile <- function(type, n, alpha=0.05) {

	# statistic upper tail percentage points
	# see D'Agostino and Stephens "Goodness-of-fit techniques", 1986

if (!is.element(type, c("greenwood", "qm", "ks", "V", "cvm"))) stop("The goodness-of-fit test must be one of: Greenwood, Quesenberry-Miller, Kolmogorov-Smirnov, V = (D+) +  (D-), or Cramer-Von Mises")

if (length(alpha)!=1) {
	stop("At this stage, alpha must be a single number")
} else if (!is.element(alpha, c(0.1, 0.05, 0.025, 0.01))) {
	stop("Significance level must be a single number among 0.1, 0.05, 0.025 or 0.01")
}

index <- is.element(c(0.1, 0.05, 0.025, 0.01), alpha)

if (type=="greenwood") {
	# quantile of n*G(n)
	if (n>500) {
		q <- qnorm(1-alpha/2, mean=2*n/(n+2), sd=sqrt(4/n))
	} else {
                # greenwood.table is loaded via sysdat.rda
		f <- approxfun(greenwood.table[,1], greenwood.table[,c(2,3,4,5)[index]])
		q <- f(n)
	}
	q <- q/n		# quantile of G(n)
} else if (type=="QM") {
	stop("no default value for Quesenberry and Miller statistic")
} else if (type=="ks") {
	q <- c(1.224, 1.358, 1.480, 1.628)[index]/(sqrt(n) + 0.12 + 0.11/sqrt(n))
} else if (type=="V") {
	q <- c(1.620, 1.747, 1.862, 2.001)[index]/(sqrt(n) + 0.155 + 0.24/sqrt(n))
} else if (type=="cvm") {
	q <- c(0.347, 0.461, 0.581, 0.743)[index]/(1+1/n) + 0.4/n - 0.6/n^2
}
	
return(q)
	
}
