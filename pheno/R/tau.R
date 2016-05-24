# Kendall's normalized tau for time series x
# modified after (Kendall and Gibbons 1990)
# input: vector x of the time series
# output : tau
tau <- function(x) {
	n <- length(x)
	if(n < 2)
	    stop("tau: not enough finite observations")
	rankx <- rank(x)
	S <- 0
	for (i in 2:n-1) {
		for (j in i:n) {
			S <- S + sign(rankx[j]-rankx[i])
		}
	}
	ntg <- length(unique(rankx)) # number of tied groups
	b <- rep(1,ntg)
	j <- 0
	for (i in 1:n) { # number of values in each tied group
		if(duplicated(sort(rankx))[i]) {  b[j] <- b[j]+1 }
		else { j <- j+1 }
	}
	tmp <- 0
	for (i in 1:ntg) { tmp <- tmp + b[i]*(b[i]-1)*(2*b[i]+5) }
	# variance of S
	var <- 1/18*(n*(n-1)*(2*n+5)-tmp)
	if(S == 0) { t <- 0 }
	else {
		if(S > 0) { t <- (S-1)/sqrt(var) }
		else  	  { t <- (S+1)/sqrt(var) }
	}
	return(t)
}
