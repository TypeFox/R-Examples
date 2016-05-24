calculate.all.pairwise.Fst <-
function(allele.counts,sample.sizes){
		number.of.populations <- nrow(allele.counts)
		list.of.pairwise.comparisons <- combn(1:number.of.populations,2)
		pairwise.Fst.matrix <- matrix(0,nrow=number.of.populations,ncol=number.of.populations)
			Fst.vector <- numeric(ncol(list.of.pairwise.comparisons))
				for(i in 1:ncol(list.of.pairwise.comparisons)){
					pair.of.allele.counts <- allele.counts[list.of.pairwise.comparisons[,i],]
					pair.of.sample.sizes <- sample.sizes[list.of.pairwise.comparisons[,i],]
					Fst.vector[i] <- calculate.pairwise.Fst(pair.of.allele.counts,pair.of.sample.sizes)
				}
			pairwise.Fst.matrix[lower.tri(pairwise.Fst.matrix)] <- Fst.vector				
			pairwise.Fst.matrix[upper.tri(pairwise.Fst.matrix)] <- t(pairwise.Fst.matrix)[upper.tri(pairwise.Fst.matrix)]
		return(pairwise.Fst.matrix)
	}
