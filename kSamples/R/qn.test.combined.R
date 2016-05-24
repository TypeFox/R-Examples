qn.test.combined <-
function (...,data = NULL, test = c("KW","vdW","NS"),
	method=c("asymptotic","simulated","exact"),
	dist=FALSE,Nsim=10000) 
{
#############################################################################
# This function qn.test.combined combines several QN K-sample test statistics 
# QN.m, m = 1,...,M, into one overall test 
# statistic QN.combined as suggested for the Kruskal-Wallis test in 
# Lehmann, E.L. (2006), Nonparametrics, Statistical Methods Based on Ranks, 
# Ch. 6, Sec. 5D.
# See also the documentation of qn.test for the comparison of a single 
# set of K samples.
# Each application of the QN K-sample test can be based on a different K > 1. 
# This combined version tests the hypothesis that all the hypotheses 
# underlying the individual K-sample tests are 
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
#        or a data frame with first column representing the responses y,
#           the second column a factor g the levels of which are used to 
#           indicate the samples within each block, and the third column 
#           a factor b indicating the block
#
#        or a formula y ~ g | b with y, g, b as in the previous situation
#
#        or just the three vectors y, g, b in this order with same meaning.
#
# 		test:	specifies the ranks scores to be used, averaging the scores
#               of tied observations. 
#               test = "KW" uses scores 1:N ( ==> Kruskal-Wallis test)
#               test = "vdW" uses the van der Waerden scores qnorm(1:N/(N+1))
#               test = "NS" uses normal scores, expected standard normal order 
#               statistics, uses function normOrder of package SuppDists.
#               Other scores could easily be added to this function and
#               to qn.test.
#
#		method 
#       can take one of three values "asymptotic", "simulated",
#		and "exact", which determines the mode of P-value calculation.
#       The asymptotic P-value based on the chi-square approximation
#       is always returned.
#       The simulated P-value simulates splits of the pooled samples in 
#		the i-th list L.i into K.i samples of sizes n.i[1], ..., n.i[K.i],
#       computing the corresponding QN statistic QN.i,
#       doing this independently for i = 1,...,M and adding the QN.i
#       to get the combined statistic QN.comb. This is repeated Nsim 
#       times and the simulated P-value is the proportion of these 
#       values that are >= the observed combined value.
#       The exact P-value should only be attempted for small M and small
#       sample sizes and requires that Nsim be set to >= the total number
#       of QN.comb enumerations. Otherwise Nsim simulations are run
#       to get a simulated P-value, as described above.
#       As example consider: M=2 with K.1 = 2, n.1[1] = 5, n.1[2] = 6,
#		K.2 = 2, n.2[1] = 5, n.2[2] = 7, then we would have
#       choose(11,5) = 462 splits of the first set of pooled samples
#       and choose(12,5) = 792 splits of the second set of pooled samples
#       and thus 462 * 792 = 365904 combinations of QN.1+QN.2 = QN.comb.
#       Thus we would need to set Nsim >= 365904 to enable exact 
#       exact enumeration evaluation of the P-value. Since these enumerated
# 		values of QN.comb need to be held inside R in a single vector,
#  		we are limited by the object size in R. In simulations the length
#       of the simulated vector of QN.comb is only Nsim and is manageable.
#
#       dist
#       takes values FALSE (default) or TRUE, where TRUE enables the 
#	return of the simulated or exact vectors of generated values
#       of QN.comb.
#       Otherwise NULL is returned for both versions
#
#       Nsim 
#       = 10000 (default), number of simulations as discussed above.
#       
#       
#
# An example:
# x1 <- c(1, 3, 2, 5, 7)
# x2 <- c(2, 8, 1, 6, 9, 4) 
# x3 <- c(12,  5,  7,  9, 11)
# y1 <- c(51, 43, 31, 53, 21, 75) 
# y2 <- c(23, 45, 61, 17, 60)
# then 
# set.seed(2627)
# qn.test.combined(list(x1,x2,x3),list(y1,y2),test="KW", method="simulated",Nsim=100000) 
# or equivalently qn.test.combined(list(list(x1,x2,x3),list(y1,y2)),
#	test="KW",method="simulated",Nsim=100000)
# produces the outcome below.
##########################################################################
# 
# Combination of Kruskal-Wallis K-Sample Tests.
# 
# Number of data sets = 2 
# 
# Sample sizes within each data set:
# Data set 1 :  5 6 5
# Data set 2 :  6 5
# Total sample size per data set: 16 11 
# Number of unique values per data set: 11 11 
# 
# Null Hypothesis:
# All samples within a data set come from a common distribution.
# The common distribution may change between data sets.
# 
# Based on Nsim = 1e+05 simulations
# for data set 1 we get
#               QN  asympt. P-value     sim. P-Value 
#       5.64851852       0.05935261       0.05156000 
# 
# for data set 2 we get
#               QN  asympt. P-value     sim. P-Value 
#       0.03333333       0.85513214       0.93001000 
# 
# Combined Kruskal-Wallis Criterion: QN.combined = QN.1+QN.2 
# 
# Based on Nsim = 1e+05 simulations
#          QN.comb  asympt. P-value     sim. P-Value 
#        5.6818519        0.1281575        0.1216000 
#
###############################################################################
# For out.qn.combined <- qn.test.combined(list(x1,x2,x3),list(y1,y2))
# or out.qn.combined <- qn.test.combined(list(list(x1,x2,x3),list(y1,y2)))
# we get the object out.qn.combined of class ksamples with the following 
# components
# > names(qn.combined.out)
# [1] "test.name" "M"         "n.samples" "nt"        "n.ties"    "qn.list"  
# [7] "qn.c"      "warning"   "null.dist" "method"    "Nsim" 
# where 
# test.name 	"Kruskal-Wallis", "van der Waerden scores", or "normal scores"
# M 			number of sets of samples being compared
# n.samples		= is a list of the vectors giving the sample sizes for each 
#             	set of samples being compared
# nt  			vector of total sample sizes involved in each of the M comparisons
# n.ties  		vector of lenth M giving the number of ties in each comparison group
# qn.list 	 	list of M data frames for the results for each of the test results
#           	corresponding to the M block
# qn.c  		2 (or 3) vector containing QN.obs, asymptotic P-value, 
#          		(and simulated or exact P-value) for the combined test.
#         		Here qn.obs is the observed value of QN.comb and 
#				P-value = P(QN.comb >= qn.obs).
# warning  	logical indicator, warning = TRUE when at least one of 
#           	the sample sizes is < 5.
# null.dist 	vector of simulated or fully enumerated QN statistics, if requested,
#               otherwise it is NULL
# method  		one of the following values: "asymptotic", "simulated", "exact"
# 				as it was ultimately used.
# Nsim  		number of simulations used, when applicable.
# Nsim
#
# Fritz Scholz, April 2012
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
# if not already in this form. It drops blocks (sublist) with at
# most one sample in it.
data.sets <- io2(...,data = data)
data.sets <- test.list(data.sets)
# end of data.sets list conversion
test <- match.arg(test)
method <- match.arg(method)
n.sizes <- NULL
M <- length(data.sets) # number of data sets

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
if(ncomb > Nsim & method == "exact") method <- "simulated"	

if( na.t > 1) cat(paste("\n",na.t," NAs were removed!\n\n"))
if( na.t == 1) cat(paste("\n",na.t," NA was removed!\n\n"))

warning <- min(n.sizes) < 5 # set warning flag for caution on
                            # trusting asymptotic p-values

# Initializing output objects
QNobs <- 0
sig <- NULL
n.ties <- NULL
nt <- NULL
mu <- NULL
qn.list <- list()
mu.c <- 0
null.dist <- NULL
if(method == "asymptotic"){dist0 <- FALSE}else{dist0 <- TRUE}
# the following loop aggregates the (estimated or exact)
# convolution distribution of the combined QN statistic versions
for(i in 1:M){
	out <- qn.test(data.sets[[i]],test=test,method=method,dist=dist0,Nsim=Nsim)
	if(dist0==T){
	    if(i == 1){
			null.dist <- out$null.dist
		}else{
			if(method=="simulated"){
				null.dist <- null.dist+out$null.dist
			}else{
				null.dist <- convvec(null.dist,out$null.dist)
			}
		}
	}
	qn.list[[i]] <- out$qn
	mu <- c(mu, length(data.sets[[i]])-1)
	QNobs <- QNobs + out$qn[1] # aggregates combined QN stats 
	n.ties <- c(n.ties, out$n.ties)
	nt <- c(nt, sum(out$ns)) # accumulates total observations in data sets
}
# get exact or simulated P-value
if(dist0==T){
    nrow <- length(null.dist)
	pval <- sum(null.dist >= QNobs)/nrow
}

# get asymptotic P-value
pval.asympt <- 1-pchisq(QNobs, sum(mu))

    if(method=="asymptotic"){
     	qn.c <- c(QNobs,pval.asympt)
	}else{
		qn.c <- c(QNobs,pval.asympt,pval)
	}
    if(method=="asymptotic"){
		names(qn.c) <- c("comb.statistic"," asympt. P-value")
    }
    if(method=="exact"){
		names(qn.c) <- c("comb.statistic"," asympt. P-value","exact P-Value")
    }
    if(method=="simulated"){
		names(qn.c) <- c("comb.statistic"," asympt. P-value","sim. P-Value")
    }
if(dist==FALSE){
	null.dist <- NULL

}
    if(test == "vdW") test.name <- "van der Waerden scores"
    if(test == "NS") test.name <- "normal scores"
    if(test == "KW") test.name <- "Kruskal-Wallis"


object <- list(test.name =test.name, M = M,
				n.samples=n.samples, nt=nt, n.ties=n.ties, qn.list=qn.list,
           		qn.c = qn.c, warning=warning,null.dist=null.dist,
				method=method,Nsim=Nsim)
class(object) <-  "kSamples"
object
}

