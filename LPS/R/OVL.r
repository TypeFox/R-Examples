## Overlapping coefficient (area of the overlap, to divide by 2 for proportion of both curves)
## Author : Sylvain Mareschal <maressyl@gmail.com>

### ## Inspired from : http://stats.stackexchange.com/questions/12209/percentage-of-overlapping-regions-of-two-normal-distributions
### OVL <- function(mu1, mu2, sd1, sd2, raw=FALSE) {
### 	min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
### 		f1 <- dnorm(x, mean=mu1, sd=sd1)
### 		f2 <- dnorm(x, mean=mu2, sd=sd2)
### 		pmin(f1, f2)
### 	}
### 	out <- integrate(min.f1f2, -Inf, Inf, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)
### 	if(isTRUE(raw)) return(out)
### 	else            return(out$value)
### }

OVL <- function(means, sds, cutoff=1e-4, n=1e4) {
    # Points to evaluate
	from <- min(qnorm(cutoff, mean=means, sd=sds, lower.tail=TRUE))
	to <- max(qnorm(cutoff, mean=means, sd=sds, lower.tail=FALSE))
	x <- seq(from=from, to=to, length.out=n)
	
	# Normal densities
	d1 <- dnorm(x, mean=means[1], sd=sds[1])
	d2 <- dnorm(x, mean=means[2], sd=sds[2])
	
	# Score
	out <- sum(pmin(d1, d2)) / sum(d1+d2)
	return(out)
}
