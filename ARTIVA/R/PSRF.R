PSRF <-
function(MCMCseqs,m){

# compute the potential scale reduction factor (PSRF) from multiple MCMC sequences (for a scalar or vector variable)
# INPUT:
#  - MCMCseqs : a mn x k matrix where m is the number of MCMC sequences and k is the dimension of the variable whose distribution is sampled
#    the m sequences are diplasyed one underneath the other
#  - m: the number of MCMC sequences
# OUPUT: PSRF
# depends on: .

	n=dim(MCMCseqs)[1]/m
	
	posJ= seq(1,dim(MCMCseqs)[1],by=n)
	phi_seq_bar = matrix(0, length(posJ),  dim(MCMCseqs)[2])
	s2_seq= matrix(0, length(posJ),  dim(MCMCseqs)[2])

	# compute within sequence mean phi_seq_bar & variance s2_seq
	for (i in 1:(length(posJ))){
		phi_seq_bar[i,]=apply(MCMCseqs[posJ[i]:(posJ[i]+n-1),],2,mean)
		s2_seq[i,]=apply(MCMCseqs[posJ[i]:(posJ[i]+n-1),],2,var)
	}
	
	# compute between sequence mean phi_bar
	phi_bar=apply(phi_seq_bar,2, mean)

	# compute within sequence variance W
	W=mean(apply(s2_seq,2,mean))
	
	# compute between sequence variance B
	B=sum((phi_seq_bar-t(matrix((phi_bar),dim(phi_seq_bar)[2],dim(phi_seq_bar)[1])))^2)*n/(dim(phi_seq_bar)[1]-1)/dim(phi_seq_bar)[2]

    # unbiased estimated variancevarPhi 
	varPhi=(n-1)/n*W + B/n

    # return PSRF
	return (mean(varPhi)/mean(W))
}
