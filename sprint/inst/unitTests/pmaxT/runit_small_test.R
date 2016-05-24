test.small_test <- function(){     
# Gene expression data from Golub et al. (1999)
# To reduce computation time and for illustrative purposes, we condider only
# the first 100 genes and use the default of B=10,000 permutations.
# In general, one would need a much larger number of permutations
# for microarray data.
	
	data(golub)
	smallgd<-golub[1:100,] 
	classlabel<-golub.cl
	
# Permutation unadjusted p-values and adjusted p-values 
# for maxT and minP procedures with Welch t-statistics
#	resT<-mt.maxT(smallgd,classlabel)
#	rawp<-resT$rawp[order(resT$index)]
#	teststat<-resT$teststat[order(resT$index)]
	
	
# Note that the test statistics used in the examples below are not appropriate 
# for the Golub et al. data. The sole purpose of these examples is to 
# demonstrate the use of the mt.maxT and mt.minP functions.
	
# Permutation adjusted p-values for maxT procedure with paired t-statistics
#	classlabel<-rep(c(0,1),19)
#	mt.maxT(smallgd,classlabel,test="pairt")
	
# Permutation adjusted p-values for maxT procedure with block F-statistics
	classlabel<-rep(0:18,2)
#	mt.maxT(smallgd,classlabel,test="blockf",side="upper")
	
	
	res_from_maxT <- mt.maxT(smallgd,classlabel,test="blockf",side="upper")
	res_from_pmaxT <- pmaxT(smallgd,classlabel,test="blockf",side="upper")
	
	invisible(checkEqualsNumeric(res_from_maxT[,2], res_from_pmaxT[,2]))
	
	invisible(checkEqualsNumeric(res_from_maxT[,3], res_from_pmaxT[,3]))
	
	invisible(checkEqualsNumeric(res_from_maxT[,4], res_from_pmaxT[,4]))
}
