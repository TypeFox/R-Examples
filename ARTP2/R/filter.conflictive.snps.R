
filter.conflictive.snps <- function(sum.stat, allele.info, options){
	
  msg <- paste("Removing SNPs with conflictive allele information:", date())
  if(options$print) message(msg)
	sum.info <- sum.stat$stat
	nstudy <- length(sum.info)
	
	foo <- function(x){
	  paste(sort(x), collapse = '')
	}
	
	ref.allele <- apply(allele.info[, c('RefAllele', 'EffectAllele')], 1, foo)
	names(ref.allele) <- allele.info$SNP
	
	exc.snps <- NULL
	for(k in 1:nstudy){
		ord.allele <- apply(sum.info[[k]][, c('RefAllele', 'EffectAllele')], 1, foo)
		rs <- names(ord.allele)
		id <- which(ord.allele != ref.allele[rs])
		if(length(id) > 0){
		  exc.snps <- c(exc.snps, rs[id])
		}
	}
	
	exc.snps <- unique(exc.snps)
	exc.snps

}


