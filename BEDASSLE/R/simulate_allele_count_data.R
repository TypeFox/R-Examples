simulate_allele_count_data <-
function(allele.frequencies,sample.sizes,phi.parameter=NULL,invariant.loci.tolerance=10,allele.frequency.numerical.shift=1e-10){
		populations <- nrow(allele.frequencies)
		loci <- ncol(allele.frequencies)
		simulated.allele.counts <- matrix(0,nrow=populations,ncol=loci)

		if(is.null(phi.parameter)){
			invariant.loci <- apply(simulated.allele.counts,2,identify_invariant_loci)
			i <- 0
				while(any(invariant.loci=="TRUE") && i < invariant.loci.tolerance){
					simulated.allele.counts <- matrix(	rbinom(n=populations*loci,
														size=sample.sizes,
														prob=allele.frequencies),
												nrow=populations,ncol=loci)
					invariant.loci <- apply(simulated.allele.counts,2,identify_invariant_loci)
					i <- i + 1
				}										
		}						
		if(!is.null(phi.parameter)){
			invariant.loci <- apply(simulated.allele.counts,2,identify_invariant_loci)
			i <- 0
				while(any(invariant.loci=="TRUE") && i < invariant.loci.tolerance){
					allele.frequencies[which(allele.frequencies == 0)] <- allele.frequency.numerical.shift
					allele.frequencies[which(allele.frequencies == 1)] <- 1-allele.frequency.numerical.shift
					
					simulated.allele.counts <- matrix(	rbetabinom(n=populations*loci,
														size=sample.sizes,
														prob=allele.frequencies,
														shape1=phi.parameter*allele.frequencies,
														shape2=phi.parameter*(1-allele.frequencies)),
												nrow=populations,ncol=loci)
					invariant.loci <- apply(simulated.allele.counts,2,identify_invariant_loci)
					i <- i + 1
				}										
		}		
		if(i == invariant.loci.tolerance){
			warnings("your data matrix contains invariant loci")
		}		
		return(simulated.allele.counts)
	}
