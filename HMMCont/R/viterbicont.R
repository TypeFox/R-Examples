viterbicont <-
function(hmm)
{
	NRv<-nrow(hmm$Parameters)
	nr<-length(hmm$Observations)
	Deltas<-matrix(NA, nrow=nr, ncol=2)
	Deltas[1,1]<-(hmm$Parameters[NRv,1])*hmm$B[1,1]
	Deltas[1,2]<-(hmm$Parameters[NRv,2])*hmm$B[1,2]
	for(i in 2: (nr))
	{
		Deltas[i,1]<-(max(Deltas[(i-1),1]*(hmm$Parameters[NRv,3]), Deltas[(i-1),2]*(hmm$Parameters[NRv,5]))) * hmm$B[i,1]
		Deltas[i,2]<-(max(Deltas[(i-1),1]*(hmm$Parameters[NRv,4]), Deltas[(i-1),2]*(hmm$Parameters[NRv,6]))) * hmm$B[i,2]
	}
	Psi<-matrix(NA, nrow=nr, ncol=2)
	Psi[1,]<-0
	for(i in 2: (nr))
	{
		Psi[i,1]<-ifelse(Deltas[(i-1),1]*(hmm$Parameters[NRv,3]) > Deltas[(i-1),2]*(hmm$Parameters[NRv,5]), 1, 2)
		Psi[i,2]<-ifelse(Deltas[(i-1),1]*(hmm$Parameters[NRv,4]) > Deltas[(i-1),2]*(hmm$Parameters[NRv,6]), 1, 2)
	}
	#Q<-matrix(NA, nrow=nr, ncol=1)
	hmm$Viterbi[nr,1]<-ifelse(Deltas[nr,1] > Deltas[nr,2], 1, 2)
	
	for(i in ((nr)-1):1)
	{
		hmm$Viterbi[i,1]<-ifelse(hmm$Viterbi[i+1,1] == 1, Psi[i+1,1], Psi[i+1,2])
	}
	return(hmm)
}
