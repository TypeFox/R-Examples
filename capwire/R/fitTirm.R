fitTirm <-
function(data, max.pop){
	
	data <- classToInd(data)
	
	data <- data[which(data[,2] != 0), ]
	
	counts <- data[,2]
	
	s <- sample.size <- sum(counts) # total observations (sample)
	
	t <- sampled.ind <- length(counts) # total number of individuals
	
	cap.ind <- s/t # mean capture rate
	
	model <- "Two.innate.rates"

	
		
## likelihood function for TIRM
	
	likTirm <- function(na, nb, alpha, cia, cib, ta, tb){
		
		
		W <- sum( sapply (1:na, function(w) log(w)))
		
		X <- sum( sapply (1:(na-ta), function(x) log(x)))		##This quantity is never used as Na.mle=Ta
		
		Y <- log(alpha / ( alpha * na + nb )) * sum(cia)
		
		Z <- sum( sapply (1:nb, function(z) log(z)))
		
		U <- sum( sapply (1:(nb-tb), function(u) log(u)))
		
		V <- log(1 / ( alpha * na + nb )) * sum(cib)
		
		if (nb == tb){										##Need this condition because Nb stating value is Tb
			
			return( W + Y + Z + V )
			
		} else {
			
			return ( W + Y + Z - U + V )
			
		}
		
	}
	
## function to get likelihood constant
	
	getLikConTirm <- function(cia, cib, ta, tb, sample.size){
		c <- - sum(sapply(c(1:ta), function(x) log(x))) - sum(sapply(c(1:tb), function(x) log(x)))

		return(c)
		
	}
	
	
	
	if (all(unique(counts) == 1)){
		warning("Sorry, your data contains only singletons and is not informative")
		ml.pop.size <- max.pop
		
		likelihood <- likTirm(nb=(max.pop - 1), na=1, alpha=1, cia=counts[1], cib=counts[-1],
							  ta=1, tb=(length(counts) - 1))
		+ getLikConTirm(cia=counts[1], cib=counts[-1], ta=1, tb=(length(counts) - 1), sample.size=sample.size)
		
		ml.nb <- max.pop-1
		ml.na <- 1
		alpha <- 1
		
		return(list(model=model, likelihood=likelihood, ml.pop.size=ml.pop.size, ml.na=ml.na, ml.nb=ml.nb, alpha=alpha, cap.ind=cap.ind, sampled.ind=sampled.ind, sample.size=sample.size, max.pop=max.pop))
	}
	

# 1: Initialize Alpha (Note Steps Refer to Those Outlined in Appendix of Miller et al. 2005)
			above.avg <- counts[which(counts > cap.ind)] 
			mean.above <- mean(above.avg)
			below.avg <- counts[which(counts <= cap.ind)]
			mean.below <- mean(below.avg)
			
			alpha.in <- mean.above / mean.below
			
# 2: Expected Capture Counts
			nn <- t
			na <- round(nn/2)
			nb <- round(nn/2)
			alpha <- alpha.in
			
			exp.a <- s * (alpha / (alpha * na + nb))
			exp.b <- s * (1 / (alpha * na + nb))
			
# 3: Assign Capture Classes
			aa <- vector()
			bb <- vector()
			
			for (i in 1:length(counts)){
				if (counts[[i]] == 1){ # assign to class b
					bb <- c(bb, counts[[i]])				
				} else {
					tmpa <- abs(counts[[i]] - exp.a)
					tmpb <- abs(counts[[i]] - exp.b)
				
					if (tmpa < tmpb){ # assign to class a
						aa <- c(aa, counts[[i]])
					}
			
					if (tmpa >= tmpb){ # assign to class b
						bb <- c(bb, counts[[i]])
					}
				}
			}

# 4: MLE Estimation of N

	ta <- length(aa)
	na <- ta
	tb <- length(bb)
	nb <- tb
	vv <- c(nb:(max.pop - na)) # possible values of nb
	lik <- sapply(vv, function(b) likTirm(na=na, b, alpha=alpha, cia=aa, cib=bb, ta=ta, tb=tb))
	
	nb <- vv[which(lik == max(lik))]
	

# 5: Bias Adjustment of alpha using equation 3 from Miller et al. 2005

	alpha <- (nb * sum(aa)) / (na * s - na * sum(aa))
			
# 6: Repeat MLE Estimation with bias corrected alpha		
			
	vv <- c(nb:(max.pop - na)) # possible values of nb
	lik <- sapply(vv, function(b) likTirm(na=na, b, alpha=alpha, cia=aa, cib=bb, ta=ta, tb=tb))
	
	nb <- vv[which(lik == max(lik))]
		
# 7: Repeat capture class assignment
	exp.a <- s * (alpha / (alpha * na + nb))
	exp.b <- s * (1 / (alpha * na + nb))

		  aa2 <- vector() # new capture assignments
			bb2 <- vector() 
			
			for (i in 1:length(counts)){
				if (counts[[i]] == 1){ # assign to class b
					bb2 <- c(bb2, counts[[i]])				
				} else {
					tmpa <- abs(counts[[i]] - exp.a)
					tmpb <- abs(counts[[i]] - exp.b)
				
					if (tmpa < tmpb){ # assign to class a
						aa2 <- c(aa2, counts[[i]])
					}
			
					if (tmpa >= tmpb){ # assign to class b
						bb2 <- c(bb2, counts[[i]])
					}
				}
			}
			
			ta2 <- length(aa2)
			tb2 <- length(bb2)

	i <- 1
# 8: Check for convergence (have the capture classes changed)
# If not, repeat steps 4-7
	while(1){
		if (tb == tb2){
			break()
		}
		if (i == 20){
			break()
		}	
		aa <- aa2
		bb <- bb2
		
    	ta <- length(aa)
    	na <- ta
		tb <- length(bb)
		nb <- max(c(nb, tb))
		vv <- c(nb:(max.pop - na)) # possible values of nb
		lik <- sapply(vv, function(b) likTirm(na=na, b, alpha=alpha, cia=aa, cib=bb, ta=ta, tb=tb))
	
		nb <- vv[which(lik == max(lik))]
	

		alpha <- (nb * sum(aa)) / (na * s - na * sum(aa))
			
		vv <- c(nb:(max.pop - na)) # possible values of nb
		lik <- sapply(vv, function(b) likTirm(na=na, b, alpha=alpha, cia=aa, cib=bb, ta=ta, tb=tb))
	
		nb <- vv[which(lik == max(lik))]
		
		exp.a <- s * (alpha / (alpha * na + nb))
		exp.b <- s * (1 / (alpha * na + nb))

	   aa2 <- vector() # new capture assignments
	    bb2 <- vector() 
			
			for (i in 1:length(counts)){
				if (counts[[i]] == 1){ # assign to class b
					bb2 <- c(bb2, counts[[i]])				
				} else {
					tmpa <- abs(counts[[i]] - exp.a)
					tmpb <- abs(counts[[i]] - exp.b)
				
					if (tmpa < tmpb){ # assign to class a
						aa2 <- c(aa2, counts[[i]])
					}
			
					if (tmpa >= tmpb){ # assign to class b
						bb2 <- c(bb2, counts[[i]])
					}
				}
			}
			
			ta2 <- length(aa2)
			tb2 <- length(bb2)
		i <- i + 1
	}
	
	likelihood <- max(lik) + getLikConTirm(aa2, bb2, ta2, tb2, sample.size)
		
	ml.pop.size <- nb + na
	
	ml.nb <- nb
	
	ml.na <- na
	
	alpha <- alpha
		
	return(list(model=model, likelihood=likelihood, ml.pop.size=ml.pop.size, ml.na=ml.na, ml.nb=ml.nb, alpha=alpha, cap.ind=cap.ind, sampled.ind=sampled.ind, sample.size=sample.size, max.pop=max.pop))
	
}
