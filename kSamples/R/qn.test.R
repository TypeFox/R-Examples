qn.test <-
function (..., data = NULL,
		test = c("KW","vdW","NS"),
		method=c("asymptotic","simulated","exact"),
		dist=FALSE,Nsim=10000) 
{
#############################################################################
# This function "qn.test" tests whether k samples (k>1) come from a common
# continuous distribution, using the QN rank test. See Lehmann (2006),
# Nonparametrics, Statistical Methods Based on Ranks, Appendix Corollary 10.
# Ties are handled by using average rank scores.
# While the asymptotic P-value is always returned, there is the option 
# to get an estimate based on Nsim simulations or an exact value based 
# on the full enumeration distribution, provided method = "exact" is chosen
# and the number of full enumerations is <= the Nsim specified.
# If the latter is not the case, simulation is used with the indicated Nsim.
# These simulated or exact P-values are appropriate under the continuity 
# assumption or, when ties are present, they are still appropriate
# conditionally on the tied rank pattern, provided randomization took 
# place in allocating subjects to the respective samples, i.e., also
# under random sampling from a common discrete parent population.
#
#
# 
# Inputs:
#       	...:	can either be a sequence of k (>1) sample vectors,
#
#				or a list of k (>1) sample vectors,
#
#                               or y, g, where y contains the concatenated
#				samples and g is a factor which by its levels
#				identifies the samples in y,
#
#                               or a formula y ~ g with y and g as in previous case.
#
#           test:	specifies the ranks scores to be used, averaging the scores
#                   of tied observations. 
#                   test = "KW" uses scores 1:N ( ==> Kruskal-Wallis test)
#                   test = "vdW" uses the van der Waerden scores qnorm(1:N/(N+1))
#                   test = "NS" uses normal scores, expected standard normal order 
#                   statistics, uses function normOrder of package SuppDists.
#                   Other scores could easily be added to this function.
#				
#			method: takes values "asymptotic", "simulated", or "exact".
#					The value "asymptotic" causes calculation of P-values
#            		using the asymptotic chi-square approximation, always done.
#
#					The value "simulated" causes estimation of P-values
#					by randomly splitting the the pooled data into
#					samples of sizes ns[1], ..., ns[k], where
#  					ns[i] is the size of the i-th sample vector,
#					and n = ns[1] + ... + ns[k] is the pooled sample size.
#					For each such random split the QN statistic is 
#					computed. This is repeated Nsim times and the proportions
#					of simulated values >= the actually observed QN value
#					is reported as P-value estimate.
#
#                   The value "exact" enumerates all n!/(ns[1]! * ... * ns[k])
#                   splits of the pooled sample and computes the QN statistic.
#					The proportion of all enumerated QN statistics
# 					that are >= the actually observed QN value
#					is reported as exact (conditional) P-value.
#
#			dist: 	= FALSE (default) or TRUE, TRUE causes the simulated
#					or fully enumerated vector of the QN statstic to be returned
#     				as null.dist.
#
#			Nsim: 	number of simulations to perform, 
#					for method = "exact" to take hold, it needs to be at least
#					equal the number of all possible splits of the pooled
#					data into samples of sizes ns[1], ..., ns[k], where
#  					ns[i] is the size of the i-th sample vector.
#
# When there are NA's among the sample values they are removed,
# with a warning message indicating the number of NA's.
# It is up to the user to judge whether such removals make sense.
#
# An example:
# z1 <- c(0.824, 0.216, 0.538, 0.685)
# z2 <- c(0.448, 0.348, 0.443, 0.722)
# z3 <- c(0.403, 0.268, 0.440, 0.087)
# qn.test(z1,z2,z3,method="exact",dist=T,Nsim=100000)
# or 
# qn.test(list(z1,z2,z3),test="KW",method="exact",dist=T,Nsim=100000)
# which produces the output below.
#############################################################################
#
# Kruskal-Wallis  k-sample test.
# 
# Number of samples:  3
# Sample sizes: 4 4 4
# Total number of values: 12
# Number of unique values: 12
# 
# Null Hypothesis: All samples come from a common population.
# 
#               QN  asympt. P-value    exact P-Value 
#        3.5769231        0.1672172        0.1729870 
# 
# 
# Warning: At least one sample size is less than 5.
# asymptotic p-values may not be very accurate.
#
#############################################################################
# In order to get the output list, call 
# qn.out <- qn.test(list(z1,z2,z3),test="KW",method="exact",dist=T,Nsim=100000)
# then qn.out is of class ksamples and has components 
# > names(qn.out)
# [1] "test.name" "k"         "ns"        "N"         "n.ties"    "qn"       
# [7] "warning"   "null.dist" "method"    "Nsim"    
#
# where
# test.name = "Kruskal-Wallis", "van der Waerden", or "normal scores"
# k = number of samples being compared
# ns = vector of the k sample sizes ns[1],...,ns[k]
# N = ns[1] + ... + ns[k] total sample size
# n.ties = number of ties in the combined set of all n observations
# qn =  2 (or 3) vector containing the QN statistics, its asymptotic P-value,
#      	(and its exact or simulated P-value). 
# warning = logical indicator, warning = TRUE indicates that at least  
#		one of the sample sizes is < 5.    
# null.dist is a vector of simulated values of the QN statistic
# 		or the full enumeration of such values.
#		This vector is given when dist = TRUE is specified, 
# 		otherwise null.dist = NULL is returned.
# method = one of the following values: "asymptotic", "simulated", "exact"
# 			as it was ultimately used.
# Nsim = number of simulations used, when applicable.
#
# The class ksamples causes qn.out to be printed in a special output
# format when invoked simply as: > qn.out
# An example was shown above.
#
# Fritz Scholz, August 2012
#
#################################################################################

ave.score <- function(z, scores){
# This function takes a data vector z and a vector scores
# of same length and returns a vector av.scores of scores 
# using average scores for each group of tied 
# observations in z. av.scores and scores have same length.
   N <- length(z)
   rz <- rank(z)
   r.rz <- rank(z,ties.method="random")
   rz.u <- unique(rz)
   av.scores <- rep(0,N)
   for(rz.ui in rz.u){
      av.scores[rz==rz.ui] <- mean(scores[r.rz[rz
                                      ==rz.ui]])
   }
  av.scores
}
	samples <- io(...,data = data)
	test <- match.arg(test)
	method <- match.arg(method)
	out <- na.remove(samples)
	na.t <- out$na.total
	if( na.t > 1) print(paste("\n",na.t," NAs were removed!\n\n"))
	if( na.t == 1) print(paste("\n",na.t," NA was removed!\n\n"))
	samples <- out$x.new
	k <- length(samples)
	if (k < 2) stop("Must have at least two samples.")
	ns <- sapply(samples, length)
	n <- sum(ns)
	if (any(ns == 0)) stop("One or more samples have no observations.")
	x <- unlist(samples)
	if(test == "KW"){ scores.vec <- 1:n }
	if (test == "NS") {
		# if (!exists("normOrder")) library(SuppDists)
		scores.vec <- normOrder(n)
	}
	if(test == "vdW") {
		scores.vec <- qnorm((1:n)/(n + 1))
	}  
	QNobs <- 0
	pval <- 0
	rx <- ave.score(x,scores.vec)
	svar <- var(rx)
	smean <- mean(rx)
	L <- length(unique(rx))
	if(dist == TRUE) Nsim <- min(Nsim,1e8) 
	ncomb <- 1
	if( method == "exact"){
		np <- n
    		for(i in 1:(k-1)){
			ncomb <- ncomb * choose(np,ns[i])
        		np <- np-ns[i]
		}
	# it is possible that ncomb overflows to Inf
		if(!(ncomb < Inf)) stop('ncomb = Inf, method = "exact" not possible\n')
	}

	if( method == "exact" & Nsim < ncomb) {
		method <- "simulated"
	}
	if( method == "exact" & dist == TRUE ) nrow <- ncomb
   	if( method == "simulated" & dist == TRUE ) nrow <- Nsim

	if( method == "simulated" ) ncomb <- 1 # don't need ncomb anymore
	if(method == "asymptotic"){
		Nsim <- 1
		dist <- FALSE
	}
	useExact <- FALSE
	if(method == "exact") useExact <- TRUE
	if(dist==T){
		QNvec <- numeric(nrow)
	}else{
    		QNvec <- 0
	}	
	out <- .C("QNtest", pval=as.double(pval),
		Nsim=as.integer(Nsim), k=as.integer(k), 
		rx=as.double(rx), ns=as.integer(ns), 
            	useExact=as.integer(useExact),
		getQNdist=as.integer(dist),
		ncomb=as.double(ncomb),QNobs=as.double(QNobs),
		QNvec = as.double(QNvec))
	QNobs <- (out$QNobs - n*smean^2)/svar
	pval <- out$pval
	if(dist){
		QNvec <- round((out$QNvec- n*smean^2)/svar,8)
	}
	pval.asympt <- 1-pchisq(QNobs,k-1)
	if(method=="asymptotic"){
		qn <- c(QNobs,pval.asympt)
	}else{
		qn <- c(QNobs,pval.asympt,pval)
	}
	if(method=="asymptotic"){
		names(qn) <- c("test statistic"," asympt. P-value")
	}
	if(method=="exact"){
		names(qn) <- c("test statistic"," asympt. P-value","exact P-Value")
	}
	if(method=="simulated"){
		names(qn) <- c("test statistic"," asympt. P-value","sim. P-Value")
	}
	warning <- FALSE
	if(min(ns) < 5) warning <- TRUE
	if(dist == FALSE | method == "asymptotic") QNvec <- NULL
	if(test == "vdW") test.name <- "van der Waerden scores"
	if(test == "NS") test.name <- "normal scores"
	if(test == "KW") test.name <- "Kruskal-Wallis"
	object <- list(test.name = test.name,
		k = k, ns = ns, N = n, n.ties = n - L,
		qn = qn, warning = warning, null.dist = QNvec,
		method=method, Nsim=Nsim)
    	class(object) <- "kSamples"
    	object
}

