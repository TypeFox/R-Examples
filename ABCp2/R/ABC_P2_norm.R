ABC_P2_norm <-
function(n,ObsMean,M_Lo,M_Hi,SD_Lo,SD_Hi,delta,iter){
	posterior<-c()
	discard<-c()
	Norm<-c()
	Avg<-c()
	Std<-c()
	i<-1
	j<-1
	k<-1
	l<-1
	m<-1

	while(i <= iter){
		avg<-runif(1,M_Lo,M_Hi)
		std<-runif(1,SD_Lo,SD_Hi)
		while(j<=n){
		norm<-round(rnorm(1, mean=avg, sd=std))	
		if(norm>0){
			Norm[j]<-norm
			j<-j+1}
			}
		
		P2<-runif(1,0,1)
		sire2<-rbinom(n,Norm,P2)
		meanP2<-mean(sire2)
		
	if(abs(meanP2 - ObsMean)>delta){
		discard[k]<-P2
		k<-k+1
		}else	
	if(abs(meanP2 - ObsMean)<=delta){
		posterior[i]<-P2
		Avg[l]<-avg
		Std[m]<-std
		i<-i+1
		l<-l+1
		m<-m+1
		}	
	
}	
	list(posterior = posterior, Avg = Avg, Std = Std)
}
