Likelihood_counts <-
function(counts,sample_sizes,allele.frequencies){
		LnL_counts_mat <- counts*log(allele.frequencies) + (sample_sizes-counts)*log(1-allele.frequencies)
		return(LnL_counts_mat)
	}
