library(randtoolbox)

frac <- function(x) x-floor(x)

torusR <- function(n, k)
	frac(1:n %o% sqrt(get.primes(k)))

torusR(10, 5) - torus(10, 5) 
