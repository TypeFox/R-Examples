baumwelchcont <-
function(hmm)
{
	NR<-nrow(hmm$Parameters)
	no<-length(hmm$Observations)
	hmm$B[,1]<-dnorm((hmm$Observations), mean=(hmm$Parameters[NR,7]), sd=sqrt(hmm$Parameters[NR,9]))
	hmm$B[,2]<-dnorm((hmm$Observations), mean=(hmm$Parameters[NR,8]), sd=sqrt(hmm$Parameters[NR,10]))
	Alphas<-matrix(NA, nrow=no, ncol=2)
	Alphas[1,1]<-(hmm$Parameters[NR,1])*hmm$B[1,1]
	Alphas[1,2]<-(hmm$Parameters[NR,2])*hmm$B[1,2]
	for(i in 2: (no))
	{
		Alphas[i,1]<-(Alphas[(i-1),1]*(hmm$Parameters[NR,3])+Alphas[(i-1),2]*(hmm$Parameters[NR,5]))*hmm$B[i,1]
		Alphas[i,2]<-(Alphas[(i-1),1]*(hmm$Parameters[NR,4])+Alphas[(i-1),2]*(hmm$Parameters[NR,6]))*hmm$B[i,2]
	}
	Betas<-matrix(NA, nrow=no, ncol=2)
	Betas[(no),1]<-1
	Betas[(no),2]<-1
	for(j in ((no)-1):1)
	{
		Betas[j,1]<-(Betas[(j+1),1]*hmm$B[(j+1),1]*(hmm$Parameters[NR,3])+Betas[(j+1),2]*hmm$B[(j+1),2]*(hmm$Parameters[NR,4]))
		Betas[j,2]<-(Betas[(j+1),1]*hmm$B[(j+1),1]*(hmm$Parameters[NR,5])+Betas[(j+1),2]*hmm$B[(j+1),2]*(hmm$Parameters[NR,6]))
	}
	Preeps<-matrix(NA, nrow=((no)-1), ncol=5)
	for(i in 1:((no)-1))
	{
		Preeps[i,1]<-Alphas[i,1]*(hmm$Parameters[NR,3])*hmm$B[i+1,1]*Betas[i+1,1]
	}
	for(i in 1:((no)-1))
	{
		Preeps[i,2]<-Alphas[i,1]*(hmm$Parameters[NR,4])*hmm$B[i+1,2]*Betas[i+1,2]
	}
	for(i in 1:((no)-1))
	{
		Preeps[i,3]<-Alphas[i,2]*(hmm$Parameters[NR,5])*hmm$B[i+1,1]*Betas[i+1,1]
	}
	for(i in 1:((no)-1))
	{
		Preeps[i,4]<-Alphas[i,2]*(hmm$Parameters[NR,6])*hmm$B[i+1,2]*Betas[i+1,2]
	}
	Preeps[,5]<-Preeps[,1]+Preeps[,2]+Preeps[,3]+Preeps[,4]
	Xis<-matrix(NA, nrow=((no)-1), ncol=4) # empty matrix
	Xis<-Preeps[,1:4]/Preeps[,5]
	Gammas<-matrix(NA, nrow=((no)-1), ncol=2)
	Gammas[,1]<-Xis[,1]+Xis[,2]
	Gammas[,2]<-Xis[,3]+Xis[,4]
	hmm$Results[NR,1]<-Alphas[(no),1]+Alphas[(no),2]
	hmm$Results[NR,2]<-(Betas[1,1]*(hmm$Parameters[NR,1])*hmm$B[1,1])+(Betas[1,2]*(hmm$Parameters[NR,2])*hmm$B[1,2])
	hmm$Results[NR,3]<-Preeps[1,5]
	hmm$Results[NR,4]<- -2*log(hmm$Results[NR,1])+20 
	hmm$Results[NR,5]<- -2*log(hmm$Results[NR,1])+20*log(no)
	hmm$Results[NR,6]<- -2*log(hmm$Results[NR,1])+20*log(log(no))
	newrowres<-matrix(NA, ncol=6, nrow=1)
	hmm$Results<-rbind(hmm$Results, newrowres)
	newrowpar<-matrix(NA, ncol=10, nrow=1)
	newrowpar[1]<-Gammas[1,1]
	newrowpar[2]<-Gammas[1,2]
	newrowpar[3]<-(sum(Xis[,1]))/sum(Gammas[,1])
	newrowpar[4]<-(sum(Xis[,2]))/sum(Gammas[,1])
	newrowpar[5]<-(sum(Xis[,3]))/sum(Gammas[,2])
	newrowpar[6]<-(sum(Xis[,4]))/sum(Gammas[,2])
	newrowpar[7]<-sum(Gammas[,1]*hmm$Observations[-(no)]) / sum(Gammas[,1])
	newrowpar[8]<-sum(Gammas[,2]*hmm$Observations[-(no)]) / sum(Gammas[,2])
	newrowpar[9]<-sum(Gammas[,1] * (hmm$Observations[-(no)] - newrowpar[7]) * (hmm$Observations[-(no)] - newrowpar[7]) ) / sum(Gammas[,1])
	newrowpar[10]<-sum(Gammas[,2] * (hmm$Observations[-(no)] - newrowpar[8]) * (hmm$Observations[-(no)] - newrowpar[8]) ) / sum(Gammas[,2])
	hmm$Parameters<-rbind(hmm$Parameters, newrowpar)	
	return(hmm)
}
