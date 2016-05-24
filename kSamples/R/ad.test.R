ad.test <-
function (..., data = NULL,
  method=c("asymptotic","simulated","exact"),dist=FALSE,Nsim=10000) 
{
#############################################################################
# This function "ad.test" tests whether k samples (k>1) come from a common
# continuous distribution, using the nonparametric (rank) test described in
# Scholz F.W. and Stephens M.A. (1987), K-sample Anderson-Darling Tests,
# Journal of the American Statistical Association, Vol 82, No. 399, 
# pp. 918-924. 
# This test is consistent against all alternatives. 
# Ties are handled by using midranks, and according to the above
# reference two versions of the test statistic are returned.
# They are labeled version 1 and version 2, in the order introduced
# in the above reference.
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
#				
#               data:   data frame with variables usable in formula input, default = NULL.
#
#			method: takes values "asymptotic", "simulated", or "exact".
#					The value "asymptotic" causes calculation of P-values
#            		using the asymptotic approximation, always done.
#
#					The value "simulated" causes estimation of P-values
#					by randomly splitting the the pooled data into
#					samples of sizes ns[1], ..., ns[k], where
#  					ns[i] is the size of the i-th sample vector,
#					and n = ns[1] + ... + ns[k] is the pooled sample size.
#					For each such random split the AD statistics are 
#					computed. This is repeated Nsim times and the proportions
#					of simulated values >= the respective actually 
#					observed AD values are reported as P-value estimates.
#
#                   The value "exact" enumerates all n!/(ns[1]! * ... * ns[k])
#                   splits of the pooled sample and computes the respective 
#					AD statistics. The proportion of all enumerated AD statistics
# 					which are >= the respective actually observed AD values
#					are reported as exact P-values.
#
#			dist: 	= FALSE (default) or TRUE, TRUE causes the simulated
#					or fully enumerated vectors of both AD statstics to be returned
#     					as null.dist1 and null.dist2.
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
# ad.test(z1,z2,z3,method="exact",dist=T,Nsim=100000)
# or 
# ad.test(list(z1,z2,z3),method="exact",dist=T,Nsim=100000) 
# which produces the output below.
#############################################################################
#  Anderson-Darling k-sample test.
#
# Number of samples:  3
# Sample sizes:  4, 4, 4
# Number of ties: 0
#
# Mean of  Anderson-Darling  Criterion: 2
# Standard deviation of  Anderson-Darling  Criterion: 0.88133
#
# T.AD = ( Anderson-Darling  Criterion - mean)/sigma
#
# Null Hypothesis: All samples come from a common population.
#
#                AD    T.AD  asympt. P-value  exact P-value
# version 1: 2.6367 0.72238          0.18525        0.20924
# version 2: 2.6200 0.70807          0.18819        0.21703
#
#
# Warning: At least one sample size is less than 5.
# asymptotic p-values may not be very accurate.
#
#############################################################################
# In order to get the output list, call 
# ad.out <- ad.test(z1,z2,z3,method="exact",dist=T,Nsim=100000)
# then ad.out is of class ksamples and has components 
# > names(ad.out)
#  [1] "test.name"  "k"          "ns"         "N"          "n.ties"    
#  [6] "sig"        "ad"         "warning"    "null.dist1" "null.dist2"
# [11] "method"     "Nsim"
#
# where
# test.name = "Anderson-Darling"
# k = number of samples being compared
# ns = vector of the k sample sizes ns[1],...,ns[k]
# N = ns[1] + ... + ns[k] total sample size
# n.ties = number of ties in the combined set of all n observations
# sig = standard deviation of the AD statistic (for continuous population case)
# ad =  2 x 3 (or 2 x 4) matrix containing the AD statistics, 
#		standardized AD statistics, its asymptotic P-value, 
#      	(and its exact or simulated P-value), for version 1 in the first row 
#		and for version 2 in the second row.
# warning = logical indicator, warning = TRUE indicates that at least  
#		one of the sample sizes is < 5.
# null.dist1 is a vector of simulated values of the AD statistic (version 1)
# 		or the full enumeration of such values.
#		This vector is given when dist = TRUE is specified, 
# 		otherwise null.dist1 = NULL is returned.
# null.dist2 is the corresponding vector for the 2nd AD statistic version.
# method = one of the following values: "asymptotic", "simulated", "exact"
# 			as it was ultimately used.
# Nsim = number of simulations used, when applicable.
#
# The class ksamples causes ad.out to be printed in a special output
# format when invoked simply as: > ad.out
# An example was shown above.
#
# Fritz Scholz, August 2012
#
#################################################################################
	samples <- io(...,data = data)
	method <- match.arg(method)
	out <- na.remove(samples)
	na.t <- out$na.total
	if( na.t > 1) print(paste("\n",na.t," NAs were removed!\n\n"))
	if( na.t == 1) print(paste("\n",na.t," NA was removed!\n\n"))
	samples <- out$x.new
	k <- length(samples)
	if (k < 2) stop("Must have at least two samples.")
	ns <- sapply(samples, length)
	if (any(ns == 0)) stop("One or more samples have no observations.")
	x <- unlist(samples)
	n <- length(x)
	Z.star <- sort(unique(x))
	L <- length(Z.star)
	if(dist == TRUE) Nsim <- min(Nsim,1e8) 
	# limits the size of null.dist1 and null.dist2
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
	dist1 <- NULL
	dist2 <- NULL
	pv <- c(NA,NA)
	getA2mat <- dist
	useExact <- FALSE
	if(method == "exact") useExact <- TRUE
	if(getA2mat){
		a2mat <- matrix(0,nrow=nrow,ncol=2)}else{
		a2mat <- 0
	}
	ans <- numeric(2)
        pval <- numeric(2)
	out0 <- .C("adkTestStat0",ans=as.double(ans),k=as.integer(k),x=as.double(x),
		ns=as.integer(ns),Z.star=as.double(Z.star),L=as.integer(L))
	if(method != "asymptotic"){	
		out1 <- .C("adkPVal0",pval=as.double(pval), Nsim=as.integer(Nsim),k=as.integer(k),
				x=as.double(x),ns=as.integer(ns),
				zstar=as.double(Z.star),L=as.integer(L),
				useExact=as.integer(useExact),getA2mat=as.integer(getA2mat),
				ncomb=as.double(ncomb),a2mat=as.double(a2mat))
	    	pv <- out1$pval
  		if(getA2mat){
    			a2mat <- matrix(out1$a2mat, nrow=nrow, ncol=2, byrow=FALSE, 
			 		dimnames=list(NULL, c("AkN2", "AakNk2")))
    			dist1 <- round(a2mat[,1],8)
    			dist2 <- round(a2mat[,2],8)
		}
	}
	AkN2 <- out0[[1]][1]
	AakN2 <- out0[[1]][2]
	coef.d <- 0
	coef.c <- 0
	coef.b <- 0
	coef.a <- 0
	H <- sum(1/ns)
	h <- sum(1/(1:(n - 1)))
	g <- 0
	for (i in 1:(n - 2)) {
	        g <- g + (1/(n - i)) * sum(1/((i + 1):(n - 1)))
	}
	coef.a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * H
	coef.b <- (2 * g - 4) * k^2 + 8 * h * k + (2 * g - 14 * h - 
	        4) * H - 8 * h + 4 * g - 6
	coef.c <- (6 * h + 2 * g - 2) * k^2 + (4 * h - 4 * g + 6) * 
	        k + (2 * h - 6) * H + 4 * h
	coef.d <- (2 * h + 6) * k^2 - 4 * h * k
	sig2 <- (coef.a * n^3 + coef.b * n^2 + coef.c * n + coef.d)/((n - 
	        1) * (n - 2) * (n - 3))
	sig <- sqrt(sig2)
	TkN <- (AkN2 - (k - 1))/sig
	TakN <- (AakN2 - (k - 1))/sig
	pvalTkN <- ad.pval(TkN, k - 1,1)
	pvalTakN <- ad.pval(TakN, k - 1,2)
	warning <- min(ns) < 5
	if(method=="asymptotic"){
	     	ad.mat <- matrix(c(signif(AkN2,5), signif(TkN, 5), 
				signif(pvalTkN, 5), signif(AakN2,3) , 
				signif(TakN, 5), signif(pvalTakN, 5)), 
				byrow = TRUE, ncol = 3)
	}else{
  	   	ad.mat <- matrix(c(signif(AkN2,5), signif(TkN, 5), 
				signif(pvalTkN, 5),  signif(pv[1],5),signif(AakN2,3) , 
				signif(TakN, 5), signif(pvalTakN, 5), 
				signif(pv[2],5)), byrow = TRUE, ncol = 4)
	}
    	if(method=="asymptotic"){
		dimnames(ad.mat) <- list(c("version 1:","version 2:"),
						c("AD","T.AD"," asympt. P-value"))
    	}
    	if(method=="exact"){
		dimnames(ad.mat) <- list(c("version 1:","version 2:"),
						c("AD","T.AD"," asympt. P-value"," exact P-value"))
    	}
    	if(method=="simulated"){
		dimnames(ad.mat) <- list(c("version 1:","version 2:"),
						c("AD","T.AD"," asympt. P-value"," sim. P-value"))
    	}
	object <- list(test.name ="Anderson-Darling",
				k = k, ns = ns, N = n, n.ties = n - L, sig = round(sig, 5), 
				ad = ad.mat, warning = warning, null.dist1 = dist1, 
				null.dist2 = dist2, method=method, Nsim=Nsim)
	class(object) <- "kSamples"
	object
}

