calculate.pairwise.Fst <-
function(allele.counts,sample.sizes){
		raw.population.allele.frequencies <- allele.counts/sample.sizes
			missing.data.loci <- which(is.na(raw.population.allele.frequencies),arr.ind=TRUE)[,2]
		if(sum(missing.data.loci) > 0){
			allele.counts <- allele.counts[,-c(missing.data.loci)]
			sample.sizes <- sample.sizes[,-c(missing.data.loci)]
		}
		population.allele.frequencies <- allele.counts/sample.sizes
		mean.allele.frequencies <- colSums(allele.counts)/colSums(sample.sizes)
			MSP <- colSums(sample.sizes*t(apply(population.allele.frequencies,1,'-',mean.allele.frequencies)^2))
			MSG <- (1/(colSums(sample.sizes-1))) * colSums(sample.sizes*population.allele.frequencies*(1-population.allele.frequencies))
			n.c <- colSums(sample.sizes)-colSums(sample.sizes^2)/colSums(sample.sizes)
		 theta <- sum(MSP-MSG) / sum(MSP + (n.c-1)*MSG)
		 return(theta)
	}
