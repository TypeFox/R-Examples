identify_invariant_loci <-
function(allele.counts){
		there.are.invariants <- FALSE
			if(length(unique(allele.counts)) < 2){
				there.are.invariants <- TRUE
			}
		return(there.are.invariants)
	}
