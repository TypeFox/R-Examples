BB_Likelihood_counts <-
function(phi,counts,sample_sizes,allele.frequencies){
		shape1 <- phi*allele.frequencies
		shape2 <- phi*(1-allele.frequencies)
		L <- lbeta(counts+shape1,sample_sizes-counts+shape2) - lbeta(shape1,shape2)
		return(L)
	}
