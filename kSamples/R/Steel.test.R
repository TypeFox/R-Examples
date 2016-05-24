Steel.test <-
function (..., data = NULL, method=c("asymptotic","simulated","exact"),
		alternative = c("greater","less","two-sided"),
		dist=FALSE,Nsim=10000) 
{
#############################################################################
# This function "Steel.test" tests whether k-1 samples (k>1) come from the
# same continuous distribution as the control sample, taking as test statistic
# the maximum standardized Wilcoxon test statistic, or the
# minimum standardized Wilcoxon test statistic, or the
# maximum absolute standardized Wilcoxon test statistic,
# (the Steel statistic) for all two sample 
# comparisons with the control sample. See Lehmann (2006),
# Nonparametrics, Statistical Methods Based on Ranks, Chap. 5.5.
# While the asymptotic P-value is always returned, there is the option 
# to get an estimated P-value based on Nsim simulations or an exact P-value
# value based on the full enumeration distribution, provided method = "exact" 
# is chosen and the number of full enumerations is <= Nsim, as specified.
# The latter makes sure that the user is aware of the computation effort involved.
# If the number of full enumerations is > Nsim, simulation is used with the 
# indicated Nsim.
# These asymptotic, simulated or exact P-values are appropriate under the 
# continuity assumption or, when ties are present, they are still appropriate
# conditionally given the tied rank pattern, provided randomization took 
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
#		method: takes values "asymptotic", "simulated", or "exact".
#			The value "asymptotic" causes calculation of P-values based
#                       on the multivariate normal approximation for the correlated 
#                       rank sums. This is always done, but when "asymptotic" is
#                       is specified, that is all that's done. Useful for large samples,
#                  	or when a fast P-value assessment is desired.
#
#			The value "simulated" causes estimation of P-values
#			by randomly splitting the pooled data into
#			samples of sizes ns[1], ..., ns[k], where
#  			ns[i] is the size of the i-th sample vector,
#			and n = ns[1] + ... + ns[k] is the pooled sample size.
#			For each such random split the Steel statistic is 
#			computed. This is repeated Nsim times and the proportions
#		      	of simulated values exceeding the actually observed Steel value
#			in the appropriate direction is reported as P-value estimate.
#
#                   	The value "exact" enumerates all ncomb = n!/(ns[1]! * ... * ns[k]!)
#                   	splits of the pooled sample and computes the Steel statistic.
#			The proportion of all enumerated Steel statistics exceeding
# 			the actually observed Steel statistic value in the appropriate 
#			direction is reported as exact (conditional) P-value.
#   			This is only done when ncomb <= Nsim.
#
#		alternative: takes values "greater", "less", and "two-sided".
#
#                     	For the value "greater" the maximum standardized treatment 
#			rank sum is used as test statistic, using conditional means and
# 			standard deviations given the overall tie pattern among all
#			n observations for standardization. The test rejects for 
#			large values of this maximum.
#
#                     	For the value "less" the minimum standardized treatment 
#			rank sum is used as test statistic. The test rejects 
#                       for low values of minimum.
#
#                     	For the value "two-sided" the maximum absolute standardized 
#			treatment rank sum is used as test statistic. The test rejects 
#                       for large values of this maximum.
#
#
#		dist: 	= FALSE (default) or TRUE, TRUE causes the simulated
#			or fully enumerated vector of the Steel statistic 
#			to be returned as null.dist. TRUE should be used 
#			judiciously, keeping in mind the size ncomb or Nsim
#			of the returned vector.
#
#		Nsim: 	number of simulations to perform, 
#
#
# When there are NA's among the sample values they are removed,
# with a warning message indicating the number of NA's.
# It is up to the user to judge whether such removals make sense.
#
# An example: using the data from Steel's paper on the effect of 
# birth conditions on IQ
# z1 <- c(103, 111, 136, 106, 122, 114)
# z2 <- c(119, 100,  97,  89, 112,  86)
# z3 <- c( 89, 132,  86, 114, 114, 125)
# z4 <- c( 92, 114,  86, 119, 131,  94)
#set.seed(27)
#Steel.test(z1,z2,z3,z4,method="simulated",
#   alternative="less",Nsim=100000)
# or
#Steel.test(list(z1,z2,z3,z4),method="simulated",
#   alternative="less",Nsim=100000)
#produces the output below.
#
#  Steel Mutiple Wilcoxon Test: k treatments against a common control
#
#
#Number of samples:  4
#Sample sizes:  6, 6, 6, 6
#Number of ties: 7
#
#
#Null Hypothesis: All samples come from a common population.
#
#Based on Nsim = 1e+05 simulations
#
#  test statistic  asympt. P-value     sim. P-Value 
#     -1.77126551       0.09459608       0.10474000 
#############################################################################
# In order to get the output list, call 
# st.out <- Steel.test(list(z1,z2,z3,z4),method="simulated",
#   alternative="less",Nsim=100000)
# then st.out is of class ksamples and has components 
# names(st.out)
# > names(st.out)
# [1] "test.name" "k"         "ns"        "N"         "n.ties"    "st"       
# [7] "warning"   "null.dist" "method"    "Nsim"      "mu"        "sig0"     
# [13] "sig"       "tau"       "W"  
#
# where
# test.name = "Steel"
# k = number of samples being compared, including the control sample
# ns = vector of the k sample sizes ns[1],...,ns[k]
# N = ns[1] + ... + ns[k] total sample size
# n.ties = number of ties in the combined set of all N observations
# st =  2 (or 3) vector containing the Steel statistics, its asymptotic P-value,
#      	(and its exact or simulated P-value). 
# warning = logical indicator, warning = TRUE indicates that at least  
#		one of the sample sizes is < 5.    
# null.dist is a vector of simulated values of the Steel statistic
# 		or the full enumeration of such values.
#		This vector is given when dist = TRUE is specified, 
# 		otherwise null.dist = NULL is returned.
# method = one of the following values: "asymptotic", "simulated", "exact"
# 			as it was ultimately used.
# Nsim = number of simulations used, when applicable.
# mu = the vector of means for the Mann-Whitney statistics W.XY
# sig0 = standard deviation of V.0, when W.X1Xi are viewed as n.i^2 * V.0 + V.i
#        with V.0, V.1, ..., V.(k-1) independent with means zero.
# sig = vector of standard deviations ofV.1, ..., V.(k-1)
# tau = vector of standard deviations of  W.X1X2, ..., W.X1Xk
#    all these means and standard deviations are conditional on the tie
#    pattern and are either used in the standardization of the W.X1Xi
#    or in the normal approximation. 
# W = vector of Mann-Whitney statistics W.X1X2, ..., W.X1Xk
#
# The class ksamples causes st.out to be printed in a special output
# format when invoked simply as: > st.out
# An example was shown above.
#
# Fritz Scholz, May 2012
#
#################################################################################


samples <- io(..., data = data)

method <- match.arg(method)
alternative <- match.arg(alternative)
if(alternative=="greater") alt <- 1
if(alternative=="less") alt <- -1
if(alternative=="two-sided") alt <- 0
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

# reverse the sign of x to change testing alt = -1 to alt = 1
if(alt==-1){
	x <- -x
	alt <- 1
}

Wvec <- numeric(k-1)
xst <- 0
for(i in 2:k){
xst <- xst + ns[i-1]
Wvec[i-1] <- sum( rank( c(x[1:ns[1]],x[xst+1:ns[i]]) )[ns[1]+1:ns[i]])-ns[i]*(ns[i]+1)/2
}
Steelobs <- 0
pval <- 0
rx <- rank(x)




dvec <- as.vector(table(rx))
d2 <- sum(dvec*(dvec-1))
d3 <- sum(dvec*(dvec-1)*(dvec-2))
n2 <- n*(n-1)
n3 <- n2*(n-2)
n0 <- ns[1]
ni <- ns[-1]
sig02 <- (n0/12)*(1-d3/n3)
sig0 <- sqrt(sig02)
sig2 <- (n0*ni/12)*(n0+1-3*d2/n2-(n0-2)*d3/n3)
sig <- sqrt(sig2)
tau <- sqrt(ni^2 * sig02 + sig2)
mu <- n0*ni/2
L <- length(unique(rx))
# computing total number of combination splits
if(dist == TRUE) Nsim <- min(Nsim,1e8) 
# limits the size of null.dist
# whether method = "exact" or = "simulated"
ncomb <- 1
np <- n
for(i in 1:(k-1)){
	ncomb <- ncomb * choose(np,ns[i])
      	np <- np-ns[i]
}
# it is possible that ncomb overflows to Inf
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
if( method == "exact") useExact <- TRUE
if(useExact){nrow <- ncomb}
if(dist==TRUE){
	Steelvec <- numeric(nrow)
}else{
    	Steelvec <- 0
}

out <- .C("Steeltest", pval=as.double(pval),
		Nsim=as.integer(Nsim), k=as.integer(k), 
		rx=as.double(rx), ns=as.integer(ns), 
            	useExact=as.integer(useExact),
		getSteeldist=as.integer(dist),
		ncomb=as.double(ncomb), alt=as.integer(alt),
                mu=as.double(mu), tau=as.double(tau),
		Steelobs=as.double(Steelobs),
		Steelvec = as.double(Steelvec))
Steelobs <- out$Steelobs
# changes Steelobs back to make up for earlier sign change
if(alternative == "less") Steelobs <- -Steelobs
pval <- out$pval
if(dist){
	Steelvec <- out$Steelvec
	if(alternative == "less") Steelvec <- - Steelvec
}
if(dist == F | method == "asymptotic") Steelvec <- NULL
if(alternative != "less"){
	pval.asympt <- Steelnormal(mu,sig0,sig,tau,Wvec,ni,alternative)}else{
	pval.asympt <- Steelnormal(mu,sig0,sig,tau,Wvec,ni,"greater")
}
if(method=="asymptotic"){
	st <- c(Steelobs,pval.asympt)
}else{
	st <- c(Steelobs,pval.asympt,pval)
}
if(method=="asymptotic"){
	names(st) <- c("test statistic"," asympt. P-value")
}
if(method=="exact"){
	names(st) <- c("test statistic"," asympt. P-value","exact P-Value")
}
if(method=="simulated"){
	names(st) <- c("test statistic"," asympt. P-value","sim. P-Value")
}
warning <- FALSE
if(min(ns) < 5) warning <- TRUE
if(dist == FALSE) null.dist <- NULL

object <- list(test.name = "Steel", k = k, alternative = alternative,
		ns = ns, N = n, n.ties = n - L,
		st = st, warning = warning, null.dist = Steelvec,
		method=method, Nsim=Nsim, W=Wvec, mu=mu, tau=tau, sig0=sig0, sig=sig)
    class(object) <- "kSamples"
    object

}

