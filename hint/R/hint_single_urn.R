###########################################################################
#
# Functions for calculating distributions associated with drawing distinct
#  categories from a single urn containing duplicates in q<=n of its
#  categories.
#
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
#
###########################################################################


#### Main user-callable functions that do the dispatching to user-invisible worker functions below. ####
#
# Density, distribution (probability mass function), quantile, and random generation functions:


.hydist.check.params <- function(n, a, q)
	{
	if(!a>=0 || !n>=0 || !q>=0 || !a<=n || !q<=n || !q>=0){
		stop("the following constraints must be met:\nn > 0\n0 <= a <= n\n0 <= q <= n\n", call. = FALSE)
		}
	cmin <- a-min(floor(a/2),q)
	cmax <- a
	crange <- cmin:cmax
	return(crange)
	}


dhydist <- function(n, a, q, range = NULL, log = FALSE)
	{
	# range is a vector giving intersection sizes for which the user wishes to retrieve probabilities.
	crange <- .hydist.check.params(n, a, q)
	if(is.null(range)){
		range <- crange
		}
	rn <- intersect(range, crange)
	if(length(rn)==0){
		dist <- data.frame(c=range, p=rep(0,length(range)))
		}
	dist <- .hydist.dup(n, a, q, range = rn)
	if(log){
		dist[,2] <- log(dist[,2])
		}
	return(dist)
	}


phydist <- function(n, a, q, vals, upper.tail = TRUE, log.p = FALSE)
	{
	# vals are the values of v for which we want cumulative probabilities.
	dist <- dhydist(n, a, q)
	if(log.p){
		dist[,2] <- log(dist[,2])
		}
	crange <- .hydist.check.params(n, a, q)
	rn <- intersect(vals, crange)
	if(length(rn)==0){
		pp <- NULL
		for(i in 1:length(vals)){
			if(vals[i] < min(crange)){
				if(upper.tail){
					pp[i] <- 1
				}else{
					pp[i] <- 0
					}
			}else{
				if(upper.tail){
					pp[i] <- 0
				}else{
					pp[i] <- 1
					}
				}
			}
		pval <- data.frame(c=vals, cum.p=rep(pp, length(vals)))
	}else{
		pv <- NULL
		for(i in 1:length(rn)){
			if(upper.tail){
				inds <- which(dist[,1]==rn[i]):nrow(dist)
			}else{
				inds <- 1:which(dist[,1]==rn[i])
				}
			pv[i] <- sum(dist[inds, 2])
			}
		pval <- data.frame(c=rn, cum.p=pv)
		}
	if(log.p){
		pval[,2] <- log(pval[,2])
		}
	return(pval)
	}


qhydist <- function(p, n, a, q, upper.tail = TRUE, log.p = FALSE)
	{
	# p is a probability.
	if(!p>=0 || !p<=1){
		stop("p must be between 0 and 1\n", call. = FALSE)
		}
	crange <- .hydist.check.params(n, a, q)
	vals <- crange
	dist <- dhydist(n, a, q, range=vals)
	pvals <- phydist(n, a, q, upper.tail=upper.tail, vals=vals)
	inds <- which(pvals[,2]<=p)
	pv <- sum(dist[inds, 2])
	qq <- pvals[max(inds),1]
	if(log.p){
		pv <- log(pv)
		}
	ret <- data.frame(c=qq, cum.p=pv)
	return(ret)
	}


rhydist <- function(num = 5, n, a, q)
	{
	crange <- .hydist.check.params(n, a, q)
	probs <- dhydist(n, a, q, range=crange)
	samp <- sample(crange, num, prob = probs[,2], replace = TRUE)
	return(samp)
	}


#### Worker function(s): ####


# Single urn with q duplicates.
.hydist.dup <- function(n, a, q, range)
	{
	dist <- NULL
	for(c in range){
		qs <- 0
		jmax <- q-a+c
		for(j in 0:jmax){
			if((q-a+c-j+1)<=0 || (2*c-a-j+1)<=0){
				qs <- qs
			}else{
				qs <- sum(qs, exp(lgamma(q+1) - lgamma(a-c+1) - lgamma(q-a+c-j+1) - lgamma(j+1) + lgamma(n-a+c-j+1) - lgamma(n-c+1) - lgamma(2*c-a-j+1) - lgamma(n+q+1) + lgamma(n+q-a+1) + lgamma(a+1)) )
				}
			}
		dist <- append(dist, qs)
		}
	ret <- data.frame(c=range,p=dist)
	return(ret)
	}


### Simulation ###

sim.hydist <- function(n, a, sims = 10000, Na = rep(2,n))
	{
	# Na is a vector with the numbers in each category.
	sdist <- NULL; na <- NULL
	for(i in 1:n){
		na <- append(na, rep(i,Na[i]))
		}

	for(i in 1:sims){
		cat("\r   Completed...",round((i/sims)*100,3),"%",collapse="")
		s1 <- sample(na, a, replace = FALSE)
		sdist[i] <- length(unique(s1))
		}
	cat("\n")
	return(sdist)
	}




