jt.test <-
function (..., data = NULL,
		method=c("asymptotic","simulated","exact"),
		dist=FALSE,Nsim=10000) 
{
#############################################################################
# This function "jt.test" tests whether k samples (k>1) come from a common
# continuous distribution, using the Jonckheere-Terpstra rank test. See Lehmann (2006),
# Nonparametrics, Statistical Methods Based on Ranks.
# The test rejects the null hypothesis of no effect when JT is too large,
# i.e., a positive trend in the samples (in the order given) seems indicated.
#
# Ties are handled by using midranks.
# While the asymptotic P-value is always returned, there is the option 
# to get a P-value estimate based on Nsim simulations or an exact value based 
# on the full enumeration distribution, provided method = "exact" is chosen
# and the number of full enumerations is <= the Nsim specified.
# If the latter is not the case, simulation is used with the indicated Nsim.
# These simulated or exact P-values are appropriate under the continuity 
# assumption or, when ties are present, they are still appropriate
# conditionally on the tied rank pattern, provided randomization took 
# place in allocating subjects to the respective samples, i.e., also
# under random sampling from a common discrete parent population.
# However, under ties the results are only meaningful conditionally given
# the observed tie pattern.
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
#				
#               data:   data frame with variables usable in formula input, default = NULL.
#           
#				
#			method: takes values "asymptotic", "simulated", or "exact".
#				The value "asymptotic" causes calculation of P-values
#            		        using the asymptotic normal approximation, always done.
#
#				The value "simulated" causes estimation of P-values
#				by randomly splitting the pooled data into
#				samples of sizes ns[1], ..., ns[k], where
#  				ns[i] is the size of the i-th sample vector,
#				and n = ns[1] + ... + ns[k] is the pooled sample size.
#				For each such random split the JT statistic is 
#				computed. This is repeated Nsim times and the proportions
#				of simulated values >= the actually observed JT value
#				is reported as P-value estimate.
#
#                   		The value "exact" enumerates all n!/(ns[1]! * ... * ns[k])
#                   		splits of the pooled sample and computes the JT statistic
#				for each such split. 
#				The proportion of all enumerated JT statistics
# 				that are >= the actually observed JT value
#				is reported as exact (conditional) P-value.
#
#			dist: 	= FALSE (default) or TRUE, TRUE causes the simulated
#					or fully enumerated vector of the QN statstic to be returned
#     					as null.dist. The length of this vector should not exceed 1e8.
#
#			Nsim: 	number of simulations to perform, 
#					for method = "exact" to take hold, it needs to be >=
#					the number of all possible splits of the pooled
#					data into samples of sizes ns[1], ..., ns[k], where
#  					ns[i] is the size of the i-th sample vector.
#
# When there are NA's among the sample values they are removed,
# with a warning message indicating the number of NA's.
# It is up to the user to judge whether such removals make sense.
#
# An example:
# x1 <- c(1,2)
# x2 <- c(1.5,2.1)
# x3 <- c(1.9,3.1)
# jt.test(x1,x2,x3,method="exact",Nsim=90)
# or 
# jt.test(list(x1,x2,x3),method="exact",Nsim=90)
# which produces the output below.
#############################################################################
# Jonckheere-Terpstra k-sample test.
#
# Number of samples:  3
# Sample sizes:  2, 2, 2
# Number of ties: 0
#
# Null Hypothesis: All samples come from a common population.
# Alternative: Samples indicate a positive trend.
#
#   test statistic               mu              sig  asympt. P-value 
#        9.0000000        6.0000000        2.5166115        0.1166151 
#    exact P-Value 
#       0.1666667 
# 
# 
# Warning: At least one sample size is less than 5,
#   asymptotic p-values may not be very accurate.

#############################################################################
# In order to get the output list, call 
# JT.out <- jt.test(list(x1,x2,x3),method="exact",dist=T,Nsim=100000)
# then JT.out is of class "kSamples" and has components 
# > names(JT.out)
# [1] "test.name" "k"         "ns"        "N"         "n.ties"    "JT"       
# [7] "warning"   "null.dist" "method"    "Nsim"    
#
# where
# test.name = "Jonckheere-Terpstra"
# k = number of samples being compared
# ns = vector of the k sample sizes ns[1],...,ns[k]
# N = ns[1] + ... + ns[k] total sample size
# n.ties = number of ties in the combined set of all n observations
# JT =  4 (or 5) vector containing the JT statistics, its mean and standard 
# deviation, its asymptotic P-value, (and its exact or simulated P-value). 
# warning = logical indicator, warning = TRUE indicates that at least  
#		one of the sample sizes is < 5.    
# null.dist is a vector of simulated values of the JT statistic
# 		or the full enumeration of such values.
#		This vector is given when dist = TRUE is specified, 
# 		otherwise null.dist = NULL is returned.
# method = one of the following values: "asymptotic", "simulated", "exact"
# 			as it was ultimately used.
# Nsim = number of simulations used, when applicable.
#
# The class "kSamples" causes JT.out to be printed in a special output
# format when invoked simply as: > JT.out
# An example was shown above.
#
# Fritz Scholz, August 2015
#
#################################################################################
JTmusig <- function (rx,ns) 
{
# this function computes the mean and standard deviation of the 
# JT statistic when ties are present.
	N <- length(rx)
	dvec <- as.vector(table(rx))
	X1 <- N*(N-1)*(2*N+5)
	X2 <- sum(ns*(ns-1)*(2*ns+5))
	X3 <- sum(dvec*(dvec-1)*(2*dvec+5))
	A1 <- (X1-X2-X3)/72
	A2 <- sum(ns*(ns-1)*(ns-2))*sum(dvec*(dvec-1)*(dvec-2))
	A2 <- A2/(36*N*(N-1)*(N-2))
	A3 <- sum(ns*(ns-1))*sum(dvec*(dvec-1))/(8*N*(N-1))
	sig <- sqrt(A1+A2+A3)
	nmat <- outer(ns,ns,"*")
	mu <- sum(nmat[upper.tri(nmat)])/2
	list(mu=mu,sig=sig)
}
##############################################################


	samples <- io(..., data = data)
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
	rx <- rank(unlist(samples))
	JTobs <- 0
	pval <- 0
	L <- length(unique(rx)) # to count ties
	if(dist == TRUE) Nsim <- min(Nsim,1e8) 
	# limits the size of null.dist	
	ncomb <- 1
	np <- n
    	for(i in 1:(k-1)){
		ncomb <- ncomb * choose(np,ns[i])
        	np <- np-ns[i]
	}
	# it is possible that ncomb overflows to Inf
	if( method == "exact" & Nsim < ncomb ) {
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
	if(dist == TRUE){
		JTvec <- numeric(nrow)
	}else{
    		JTvec <- 0
	}

	out <- .C("JTtest", pval=as.double(pval),
			Nsim=as.integer(Nsim), k=as.integer(k), 
			rx=as.double(rx), ns=as.integer(ns), 
            		useExact=as.integer(useExact),
			getJTdist=as.integer(dist),
			ncomb=as.double(ncomb),JTobs=as.double(JTobs),
			JTvec = as.double(JTvec))
	JTobs <- out$JTobs
	musig <- JTmusig(rx,ns)
	mu <- musig$mu
	sig <- musig$sig
	pval <- out$pval
	if(dist){
		JTvec <- round(out$JTvec,8)
		}else{
		JTvec <- NULL
	}
	pval.asympt <- 1-pnorm((JTobs - mu)/sig)
	if(method=="asymptotic"){
		JT <- c(JTobs,mu,sig,pval.asympt)
	}else{
		JT <- c(JTobs,mu,sig,pval.asympt,pval)
	}
	if(method=="asymptotic"){
		names(JT) <- c("test statistic","mu","sig"," asympt. P-value")
	}
	if(method=="exact"){
		names(JT) <- c("test statistic","mu","sig"," asympt. P-value","exact P-Value")
	}
	if(method=="simulated"){
		names(JT) <- c("test statistic","mu","sig"," asympt. P-value","sim. P-Value")
	}
	warning <- FALSE
	if(min(ns) < 5) warning <- TRUE
	if(dist == FALSE) null.dist <- NULL
	test.name <- "Jonckheere-Terpstra"

	object <- list(test.name = test.name,
		k = k, ns = ns, N = n, n.ties = n - L,
		JT = JT, warning = warning, null.dist = JTvec,
		method=method, Nsim=Nsim)
    class(object) <- "kSamples"
    object

}

