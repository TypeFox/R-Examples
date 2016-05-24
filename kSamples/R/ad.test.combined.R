ad.test.combined <-
function (..., data = NULL, method = c("asymptotic","simulated","exact"),
	dist = FALSE, Nsim = 10000) 
{
#############################################################################
# This function ad.test.combined combines several Anderson-Darling 
# K-sample test statistics AD.m, m = 1,...,M, into one overall test 
# statistic AD.combined as suggested in Section 8 of Scholz F.W. and 
# Stephens M.A. (1987), K-sample Anderson-Darling Tests,
# Journal of the American Statistical Association, Vol 82, No. 399, 
# pp. 918-924.
# See also the documentation of ad.test for the comparison of a single 
# set of K samples.
# Each application of the Anderson-Darling K-sample test can be 
# based on a different K > 1. This combined version tests the hypothesis 
# that all the hypotheses underlying the individual K-sample tests are 
# true simultaneously.
# The individual K-sample test hypothesis is that all samples from 
# the m-th group come from a common population. However, that population 
# may be different from one individual K-sample situation to the next. 
# Such a combined test is useful in 
#
# 1) examining intra-laboratory measurement equivalence, when samples from 
#    the same material or batch are compared for M laboratories and such
#    comparisons are made for samples from several different materials or
#    batches and one assessment over all materials/batches is desired.
#
# 2) analyzing treatment effects in randomized complete or incomplete 
#    block designs.
#
# When there are NA's among the sample values they are removed,
# with a warning message indicating the number of NA's.
# It is up to the user to judge whether such removals make sense.
#
# Input: ...  
#        can take the form of several lists, 
#        say L.1,...,L.M, where list L.i contains 
#        K.i sample vectors of respective sizes n.i[1], ..., n.i[K.i] 
#        (n.i[j] > 4 is recommended)
#
#        or a single list of such lists
#
#        or a data frame with 3 columns, the first column representing 
#	 the responses y, the second column a factor g the levels of 
#  	 which are used to indicate the samples within each block, and 
#	 the third column a factor b indicating the block
#
#        or a formula y ~ g | b with y, g, b as in the previous situation,
#        where y, g, b may be variable names in a provided data frame, say dat,
#	 supplied via data = dat,
#
#        or just the three vectors y, g, b in this order with same meaning.
#
# data 
#  	an optional data frame containing the variable names used in formula input, 
#       default is NULL, in which case the used variables must exist in the calling
#       environment.
#
# method 
#       can take one of three values "asymptotic", "simulated",
#		and "exact", which determines the mode of P-value calculation.
#       The asymptotic P-value is always returned.
#       The simulated P-value simulates splits of the pooled samples in 
#		the i-th list L.i into K.i samples of sizes n.i[1], ..., n.i[K.i],
#       computing the corresponding AD statistic AD.i (both versions),
#       doing this independently for i = 1,...,M and adding the AD.i
#       to get the combined statistic AD.comb. This is repeated Nsim 
#       times and the simulated P-value is the proportion of these 
#       values that are >= the observed combined value.
#       The exact P-value should only be attempted for small M and small
#       sample sizes and requires that Nsim be set to >= the total number
#       of AD.comb enumerations. Otherwise Nsim simulations are run
#       to get a simulated P-value, as described above.
#       As example consider: M=2 with K.1 = 2, n.1[1] = 5, n.1[2] = 6,
#		K.2 = 2, n.2[1] = 5, n.2[2] = 7, then we would have
#       choose(11,5) = 462 splits of the first set of pooled samples
#       and choose(12,5) = 792 splits of the second set of pooled samples
#       and thus 462 * 792 = 365904 combinations of AD.1+AD.2 = AD.comb.
#       Thus we would need to set Nsim >= 365904 to enable exact 
#       exact enumeration evaluation of the P-value. Since these enumerated
# 		values of AD.comb need to be held inside R in a single vector,
#  		we are limited by the object size in R. In simulations the length
#       of the simulated vector of AD.comb is only Nsim and is manageable.
#
# dist
#       takes values FALSE (default) or TRUE, where TRUE enables the 
#		return of the simulated or exact vectors of generated values
#       for both versions of AD.comb.
#       Otherwise NULL is returned for both versions
#
# Nsim  = 10000 (default), number of simulations as discussed above.
#       
#       
#
# An example:
# x1 <- c(1, 3, 2, 5, 7), x2 <- c(2, 8, 1, 6, 9, 4), and 
# x3 <- c(12,  5,  7,  9, 11)
# y1 <- c(51, 43, 31, 53, 21, 75) and y2 <- c(23, 45, 61, 17, 60)
# then 
# set.seed(2627)
# ad.test.combined(list(x1,x2,x3),list(y1,y2),method="simulated",Nsim=100000) 
# or equivalently ad.test.combined(list(list(x1,x2,x3),list(y1,y2)),
#	method="simulated",Nsim=100000)
# produces the outcome below.
##########################################################################
# 
# Combination of Anderson-Darling K-Sample Tests.
# 
# Number of data sets = 2 
# 
# Sample sizes within each data set:
# Data set 1 :  5 6 5
# Data set 2 :  6 5
# Total sample size per data set: 16 11 
# Number of unique values per data set: 11 11 
# 
# AD.i = Anderson-Darling Criterion for i-th data set
# Means: 2 1 
# Standard deviations: 0.92837 0.64816 
# 
# T.i = (AD.i - mean.i)/sigma.i
#  
# Null Hypothesis:
# All samples within a data set come from a common distribution.
# The common distribution may change between data sets.
# 
# Based on Nsim = 1e+05 simulations
# for data set 1 we get
#               AD   T.AD  asympt. P-value  sim. P-value
# version 1: 3.316 1.4176         0.088063       0.09868
# version 2: 3.510 1.6286         0.070278       0.09115
# 
# for data set 2 we get
#                 AD     T.AD  asympt. P-value  sim. P-value
# version 1: 0.37267 -0.96786          0.96305       0.94529
# version 2: 0.33300 -1.02930          0.98520       0.93668
# 
#
# Combined Anderson-Darling Criterion: AD.comb = AD.1+AD.2 
# Mean = 3    Standard deviation = 1.13225 
# 
# T.comb = (AD.comb - mean)/sigma
# 
# Based on Nsim = 1e+05 simulations
#            AD.comb    T.comb  asympt. P-value  sim. P-value
# version 1: 3.68867 0.6082333        0.2205062       0.23733
# version 2: 3.84300 0.7445375        0.1902867       0.21825
#
#
###############################################################################
# For out.ad.combined <- ad.test.combined(list(x1,x2,x3),list(y1,y2))
# or out.ad.combined <- ad.test.combined(list(list(x1,x2,x3),list(y1,y2)))
# we get the object out.ad.combined of class ksamples with the following 
# components
# > names(out.ad.combined)
# > names(ad.c.out)
#  [1] "test.name"  "M"          "n.samples"  "nt"         "n.ties"    
#  [6] "ad.list"    "mu"         "sig"        "ad.c"       "mu.c"      
# [11] "sig.c"      "warning"    "null.dist1" "null.dist2" "method"    
# [16] "Nsim" 
# where 
# test.name = "Anderson-Darling"
# M = number of sets of samples being compared
# n.samples = is a list of the vectors giving the sample sizes for each 
#             set of samples being compared
# nt = vector of total sample sizes involved in each of the M comparisons
# n.ties = vector of lenth M giving the number of ties in each comparison group
# ad.list = list of M data frames for the results for each of the test results
#           corresponding to the M block
#
# mu = vector of means of the M AD statistics
# sig = vector of standard deviations of the M AD statistics
# ad.c = 2 x 3 (or 2 x 4) matrix containing the AD statistics, 
#		standardized AD statistics, its asymptotic P-value, 
#      	(and its exact or simulated P-value), for version 1 in the first row 
#		and for version 2 in the second row.
# mu.c = mean of the combined AD statistic
# sig.c = standard deviation of the combined AD statistic
# warning = logical indicator, warning = TRUE when at least one of 
#           the sample sizes is < 5.
# null.dist1 enumerated or simulated values of AD statistic, version 1
# null.dist2 enumerated or simulated values of AD statistic, version 2
# method the method that was used: "asymptotic", "simulated", "exact".
# Nsim the number of simulations that were used.
#
# Fritz Scholz, August 2012
#####################################################################################
convvec <- function(x1,x2){
#----------------------------------------------------
# R routine for calling convvec in conv.c
# created 02.05.2012  Fritz Scholz
#----------------------------------------------------
n1 <- length(x1)
n2 <-length(x2)
n <- n1*n2
x <- numeric(n)
 out <- .C("convvec",  x1=as.double(x1), n1=as.integer(n1), 
					x2=as.double(x2),n2=as.integer(n2), 
					x=as.double(x),n=as.integer(n))
out$x
}


# the following converts individual data sets into a list of such,
# if not already in this form.

data.sets <- io2(...,data=data)
data.sets <- test.list(data.sets)

# end of data.sets list conversion

method <- match.arg(method)
n.sizes <- NULL
M <- length(data.sets) # number of data sets (blocks)

n.data <- sapply(data.sets, length) 
		# gets vector of number of samples in each component of data.sets
n.samples <- list() # intended to contain vectors of sample sizes for 
                    # for each respective data set.
na.t <- 0 # intended to tally total of NA cases
ncomb <- 1 # intended to hold the total number of evaluations 
           # of the full convolution distribution, to check
           # whether exact computation is reasonable. 

for(i in 1:M){
	out <- na.remove(data.sets[i])
	na.t<- na.t + out$na.total
	data.sets[i] <- out$x.new # replace data set i with the one that has
                              # NA's removed
	n.sample <- sapply(data.sets[[i]], length) 
				# contains sample sizes for i-th data set
	n.sizes <- c(n.sizes, n.sample) 
				# accumulates all sample size, warning purpose only
	if(any(n.sample==0))
		stop(paste("one or more samples in data set", i,
                       "has no observations"))
	n.samples[[i]] <- n.sample 
	N <- sum(n.sample) # total observations in i-th data set
    	k <- length(n.sample) # number of samples in i-th data set

    # updating ncomb
   	ncomb <- ncomb * choose(N,n.sample[1])
	for(j in 2:k){
		N <- N-n.sample[j-1]
		ncomb <- ncomb * choose(N,n.sample[j])
	} # end of ncomb update for data set i
}
Nsim <- min(Nsim,1e7) 

if(ncomb > Nsim & method == "exact") method <- "simulated"	

if( na.t > 1) cat(paste("\n",na.t," NAs were removed!\n\n"))
if( na.t == 1) cat(paste("\n",na.t," NA was removed!\n\n"))

warning <- min(n.sizes) < 5 # set warning flag for caution on
                            # trusting asymptotic p-values


# Initializing output objects
AD <- 0
sig <- NULL
n.ties <- NULL
nt <- NULL
mu <- NULL
ad.list <- list()
mu.c <- 0
dist1 <- NULL
dist2 <- NULL
if(method == "asymptotic"){dist0 <- FALSE}else{dist0 <- TRUE}
# the following loop aggregates the (estimated or exact)
# convolution distribution of the combined AD statistic versions
for(i in 1:M){
	out <- ad.test(data.sets[[i]],method=method,dist=dist0,Nsim=Nsim)
	if(dist0==TRUE){
	    if(i == 1){
			dist1 <- out$null.dist1
			dist2 <- out$null.dist2
		}else{
			if(method=="simulated"){
				dist1 <- dist1+out$null.dist1
				dist2 <- dist2+out$null.dist2
			}else{
				dist1 <- convvec(dist1,out$null.dist1)
				dist2 <- convvec(dist2,out$null.dist2)
			}
		}
	}
	ad.list[[i]] <- out$ad
	sig.i <- out$sig
	mu <- c(mu, length(data.sets[[i]])-1)
	AD.i <- out$ad[,1]
	sig <- c(sig, sig.i) # accumulated stand. dev.'s of AD stats
	AD <- AD+AD.i # aggregates combined AD stats (version 1 and 2)
	mu.c <- mu.c + length(data.sets[[i]]) - 1 # aggregates mean of combined AD stats
	n.ties <- c(n.ties, out$n.ties)
	nt <- c(nt, sum(out$ns)) # accumulates total observations in data sets
}
AD <- as.numeric(AD)
# get exact or simulated P-value
if(dist0==T){
    	nrow <- length(dist1)
	ad.pval1.dist <- sum(dist1 >= AD[1])/nrow
	ad.pval2.dist <- sum(dist2 >= AD[2])/nrow
}
sig.c <- sqrt(sum(sig^2)) # standard deviation of combined AD stats
tc.obs <- (AD - mu.c)/sig.c # standardized values of AD stats

# get asymptotic P-value
ad.pval1 <- ad.pval(tc.obs[1], mu.c,1)
ad.pval2 <- ad.pval(tc.obs[2], mu.c,2)

    if(method == "asymptotic"){
		ad.c <- matrix(c(signif(AD[1],7),signif(tc.obs[1],7),round(ad.pval1,7),
           signif(AD[2],7),signif(tc.obs[2],7), round(ad.pval2,7)),
           byrow=TRUE, ncol=3)
		dimnames(ad.c) <- list(c("version 1:","version 2:"),
						c("AD.comb","T.comb"," asympt. P-value"))
    }
    if(method == "exact"){
		ad.c <- matrix(c(signif(AD[1],7),signif(tc.obs[1],7),round(ad.pval1,7),
		   round(ad.pval1.dist,7),
           signif(AD[2],7),signif(tc.obs[2],7), round(ad.pval2,7),
		   round(ad.pval2.dist,7)),byrow=TRUE, ncol=4)
		dimnames(ad.c) <- list(c("version 1:","version 2:"),
						c("AD.comb","T.comb"," asympt. P-value"," exact P-value"))
    }
    if(method == "simulated"){
		ad.c <- matrix(c(signif(AD[1],7),signif(tc.obs[1],7),round(ad.pval1,7),
		   round(ad.pval1.dist,7),
           signif(AD[2],7),signif(tc.obs[2],7), round(ad.pval2,7),
		   round(ad.pval2.dist,7)),byrow=TRUE, ncol=4)

		dimnames(ad.c) <- list(c("version 1:","version 2:"),
						c("AD.comb","T.comb"," asympt. P-value"," sim. P-value"))
    }

if(dist==FALSE){
	dist1 <- NULL
	dist2 <- NULL
}


object <- list(test.name ="Anderson-Darling",
				M=M, n.samples=n.samples, nt=nt, n.ties=n.ties, ad.list=ad.list,
           		mu=mu, sig=sig, ad.c = ad.c, mu.c=mu.c,
           		sig.c=round(sig.c,5), warning=warning,null.dist1=dist1,
				null.dist2=dist2,method=method,Nsim=Nsim)
class(object) <-  "kSamples"
object
}

