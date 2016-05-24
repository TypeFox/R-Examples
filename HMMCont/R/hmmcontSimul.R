hmmcontSimul <-
function(hmm, n)
{
	Nr<-nrow(hmm$Parameters)
	Obs<-matrix(NA, nrow=n, ncol=2)
	Obs[,1]<-rnorm(n, mean=(hmm$Parameters[Nr,7]), sd=(sqrt(hmm$Parameters[Nr,9])))
	Obs[,2]<-rnorm(n, mean=(hmm$Parameters[Nr,8]), sd=(sqrt(hmm$Parameters[Nr,10])))
	States<-matrix(0, nrow=n, ncol=1)
	States[1,1]<-rbinom(1, size=1, prob=(hmm$Parameters[Nr,2]))
	for(i in 2:n)
	{
		States[i,1]<-ifelse((States[(i-1),1]==0), rbinom(1, size=1, prob=(hmm$Parameters[Nr,4])), rbinom(1, size=1, prob=(hmm$Parameters[Nr,6])))
	}
	MarkovChain<-(States+1)	
	SimulatedObservation<-matrix(NA, nrow=n, ncol=1)
	for(j in 1:n)
	{
		SimulatedObservation[j,1]<-ifelse(MarkovChain[j,1] == 1, Obs[j,1], Obs[j,2])
	}
	hmmsimul<-list(StateProcess1=(Obs[,1]), StateProcess2=(Obs[,2]), MarkovChain=(MarkovChain), SimulatedObservation=(SimulatedObservation))
	class(hmmsimul)<-"SimulContHMM"
	return(hmmsimul)
}
