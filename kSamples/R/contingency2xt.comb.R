contingency2xt.comb <-
       function (..., method=c("asymptotic","simulated","exact"),
                  dist=FALSE,Nsim=10000) 
{
#################################################################
# This function contingency2xt.comb combines several 2 x t 
# contingency table analyses over M blocks (possibly different t 
# across blocks) by adding the Kruskal-Wallis test criteria for
# each of the blocks.
# This follows the same pattern as generally combinging 
# Kruskal-Wallis tests across blocks as suggested in 
# Lehmann, E.L. (2006), Nonparametrics, Statistical Methods Based 
# on Ranks, Ch. 6, Sec. 5D.
# See also the documentation of contingency2xt for the analysis
# of a single 2 x t table of counts.
# This combined version tests the hypothesis that all the 
# hypotheses underlying the individual contingency tables are 
# true simultaneously and is relevant when randomizations or 
# samples are independent from block to block
#
# Input: ...  
#        can take the form of several lists, 
#        say L.1,...,L.M, where list L.i contains 
#        2 sample vectors of length t > 1 each, but the
#	    t may vary from list to list,	
#
#        or a single list of such lists.
#
#	method 
#       can take one of three values "asymptotic", "simulated", 
#       and "exact", 
#	which determines the mode of P-value calculation.
#       The asymptotic P-value, based on the chi-square 
#       approximation, is always returned.
#       The simulated P-value simulates counts for the tables 
#       conditioned on the observed marginal totals 
#       (see contingency2xt), doing this independently across
#       tables, and computing the corresponding 
#       Kruskal-Wallis statistics KW.i, for i = 1,...,M.
#       Adding the KW.i yields the combined statistic
#       KW.combined. 
#   	   This is repeated Nsim times and the simulated P-value is
#       the proportion of these values that are >= the observed
#       combined value.
#       The exact P-value should only be attempted for small M
#       and small marginal totals and requires that Nsim be set
#       to >= the total number of KW.combined enumerations. 
#       Otherwise Nsim simulations are run to get a simulated 
#       P-value, as described above.
#       As example consider: M=2 with t.1 = 3 columns in the
#       first table and row totals m.1 = 40, n.1[2] = 60, while
#       the second table has t.2 = 4 columns, with row totals 
#       m.2 = 30 and n.2 = 25.
#       Then we would have choose(40+2,2) = 861 possible counts
#       for table 1
#       and choose(30+3,3) = 5456 counts for table 2, thus 
#       861*5456 = 4697616 possible count configurations for both
#       tables jointly. Thus one should choose Nsim >=  4697616.
#       However, the ultimate distributions KW.combined may have 
#       far fewer unique values.
#
#     dist
#       takes values FALSE (default) or TRUE, where TRUE enables
#       the return of the simulated or exact distributions of 
#       KW.combined. Otherwise NULL is returned for both versions
#
#     Nsim 
#       = 10000 (default), number of simulations as discussed 
#       above.
#       
#       
#
# An example:
# contingency2xt.comb(list(c(15,12,25),c(12,5,7)), 
#                     list(c(12,6,4),c(6,12,3)),
#                     method="exact",dist=F, Nsim=1e6)
# produces the outcome below.
#################################################################
# 
#   Combined Kruskal-Wallis Tests for 2 x t Contingency Tables
#
# for data set 1 we get
#       observed KW   asympt. P-value     exact P-Value 
#         3.4539090         0.1778252         0.1786888 
#
# for data set 2 we get
#       observed KW   asympt. P-value     exact P-Value 
#         4.0259740         0.1335890         0.1673558 
# 
# Combined Criterion: KW.combined = KW.1+KW.2 
# 
#      KW.combined  asympt. P-value    exact P-Value 
#        7.4798830        0.1125996        0.1099804 
#
#
#################################################################
# For
# out <- contingency2xt.comb(list(c(15,12,25),c(12,5,7)), 
#                            list(c(12,6,4),c(6,12,3)),
#                            method="exact",dist=F, Nsim=1e6)
# we get the object out of class ksamples with the following 
# components
#  names(out)
# [1] "test.name"  "t"  "M" "kw.list"   "kw.c"  "null.dist"
# [7] "method"    "Nsim"  
# where 
# test.name 	= "Combined 2 x t Contingency Tables"
# t		= vector giving the number of columns for each table
# M  		= number of tables
# kw.list 	= list of M vectors holding the results for each of 
#            the tests corresponding to the M blocks
# kw.c  	= 2 (or 3) vector containing the observed KW.combined, 
#            asymptotic P-value, (and simulated or exact P-value) 
#            for the combined test.
# null.dist 	= L x 2 matrix (if dist = TRUE), with first 
#                 column holding the unique, simulated or fully 
#                 numerated KW statistics, and the second column 
#		       holding the corresponding relative frequencies 
#                 or probabilities.
#                 If dist = FALSE we get null.dist = NULL.
# method  	  one of the following values: "asymptotic",
#              "simulated", "exact" as it was ultimately used.
# Nsim  		number of simulations used, when applicable.
#
# Fritz Scholz, April 2012
#################################################################

# the following converts individual data set lists into a list of
# such, if not already in this form.
if(length(list(...)) == 1) {
	if(is.list(...) & is.list(...[[1]])){
        	data.sets <- list(...)[[1]]}else{
       		stop("you need more than 1 block of data sets\n")
	}
   	}else {
        	data.sets <- list(...)
	} 
# end of data.sets list conversion
method <- match.arg(method)
M <- length(data.sets) # number of data sets
if(M < 2) 
	stop("To combine test results you must have at least two data sets.")
tvec <- numeric(M)
ncomb <- 1
for(i in 1:M){
	n.sample <- sapply(data.sets[[i]], length) 
        if(n.sample[1] != n.sample[2])  stop("Not all count vectors have the same length.")
        tvec[i] <- n.sample[1]
	m <- sum(data.sets[[i]][[1]])
	ncomb <- ncomb * choose(m+tvec[i] - 1,m) 
}
if(ncomb > Nsim & method == "exact") method <- "simulated"	

# Initializing output objects
KWobs <- 0
null.dist <- NULL
kw.list <- list()

if(method == "asymptotic"){
	dist0 <- FALSE
	for(i in 1:M){
		out <- contingency2xt(data.sets[[i]][[1]],data.sets[[i]][[2]],
			method=method,dist=dist0,tab0=FALSE,Nsim=Nsim)
		kw.list[[i]] <- out$KW.cont
		KWobs <- KWobs + out$KW.cont[1] 
                  # aggregates combined KW stats
	}
}
# the following loops aggregate the (estimated or exact)
# convolution distribution of the combined KW statistics
if(method == "simulated"){
	dist0 <- TRUE
	for(i in 1:M){
		out <- contingency2xt(data.sets[[i]][[1]],data.sets[[i]][[2]],
			method=method,dist=dist0,tab0=FALSE,Nsim=Nsim)
		if(i == 1){
			null.dist <- out$null.dist
			}else{
				null.dist <- null.dist + out$null.dist 
			}
		kw.list[[i]] <- out$KW.cont
		KWobs <- KWobs + out$KW.cont[1] 
               # aggregates combined KW stats 
	}
}
if(method == "exact"){
	dist0 <- TRUE
	for(i in 1:M){
		out <- contingency2xt(
                   data.sets[[i]][[1]],data.sets[[i]][[2]],
			   method=method,dist=dist0,Nsim=Nsim)
		if(i == 1){
			null.dist <- out$null.dist
			}else{
				null.dist <-   
                           conv(null.dist[,1],null.dist[,2],
					out$null.dist[,1],out$null.dist[,2])
			}
		kw.list[[i]] <- out$KW.cont
		KWobs <- KWobs + out$KW.cont[1] 
               # aggregates combined KW stats 
	}
}


# get exact or simulated P-value
if(method == "simulated"){
	pval <- sum(null.dist >= KWobs)/Nsim
}
if(method == "exact"){
	pval <- sum(null.dist[null.dist[,1] >= KWobs,2])
}

# get asymptotic P-value
pval.asympt <- 1-pchisq(KWobs, sum(tvec) - M)

    if(method=="asymptotic"){
     	kw.c <- c(KWobs,pval.asympt)
	}else{
		kw.c <- c(KWobs,pval.asympt,pval)
	}
    if(method=="asymptotic"){
		names(kw.c) <- c("KW.combined"," asympt. P-value")
    }
    if(method=="exact"){
		names(kw.c) <- 
              c("KW.combined"," asympt. P-value","exact P-Value")
    }
    if(method=="simulated"){
		names(kw.c) <- 
              c("KW.combined"," asympt. P-value","sim. P-Value")
    }
if(dist==FALSE){
	null.dist <- NULL
}
if(method== "simulated" & dist==TRUE){
	out <- table(round(null.dist,6))
        null.dist <-  
           cbind(as.numeric(names(out)),as.numeric(out)/Nsim)
	  dimnames(null.dist) <- list(NULL,c("KW","prob"))
}

object <- list(test.name =paste(
            "Combined 2 x t Contingency Tables"), 
            t = tvec, M = M, kw.list = kw.list, kw.c = kw.c, 
            null.dist = null.dist, method = method, Nsim = Nsim)
class(object) <-  "kSamples"
object

}

