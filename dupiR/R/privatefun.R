##Wrote by Federico Comoglio (D-BSSE, ETH Zurich) and Maurizio Rinaldi (Dipartimento di Scienza del Farmaco, University Piemonte Orientale 'A. Avogadro') 
##Last update: 27/05/2013


################################
#Private functions, not exported
################################

F.aux <- function(k,n,x,t) {	#k=counts, n=total number of objects, x = 1-r, t = an index
	T <- cumsum(t)	#the cumulative sum of t_i
	Tm <- sum(t)	#the total T
	m <- length(t)	#the number of measurements
	term <- prod(choose(n + c(0, T)[-(m + 1)], k - t)) * prod(choose(T,t)) * x^(n + Tm) / (1 - x)^(1 + Tm)
	return(term)
}

#Function F
F <- function(k,n,x) {		#k=counts, n=total number of objects, x = 1-r
	I <- list()
	m <- length(k)  #the number of measurements
	for (i in 1 : m) I[[i]] <- 0 : k[i]	#element i of I is the set of indices for summation i
	sum_indices <- expand.grid(I) 	#takes cartesian product I_1 x ... x I_m
	tot <- sum( apply(sum_indices, 1, F.aux, k = k, n = n , x = x) )
	return (tot)
}

normalizeConstant <- function(X, k, n1, n2) {	#n1,n2 extremes of the DUP support, X = prod(x_i)
	F(k, n1, X) - F( k, n2 + 1, X)
}

#Posterior distribution with replacement
getPwithR <- function(n, k.vec, X, denominator) {
	numerator <- X ^ n * prod( choose(n, k.vec) )
	prob <- numerator / denominator
	return(prob)
}

getECDF <- function(posterior)
	return(cumsum(posterior))

Clough <- function(object, n1, n2, b = 1e-10) {
	a <- 1
	k.vec <- object@counts
	r.vec <- object@fractions
	#choose support
	if(missing(n1) & missing(n2)) {
		n1 <- object@n1
		n2 <- object@n2
	}
	if((n2-n1+1) < 1e5) {
		s <- n1 : n2
		posterior <- dgamma(s, a + sum(k.vec), b + sum(r.vec))
		return(posterior)
	}
	else {  #if n1:n2 has more than 1e4 points we use dgamma
		object@gamma <- TRUE
		return(NULL)
	}
}
