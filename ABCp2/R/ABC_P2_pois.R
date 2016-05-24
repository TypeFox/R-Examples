ABC_P2_pois <-
function(n,ObsMean, L_Lo, L_Hi,delta,iter){
	posterior<-c()
	discard<-c()
	pois<-c()
	Lambda<-c()

	i<-1
	j<-1
	k<-1
	l<-1

	while(i <= iter){
		lambda<-runif(1,L_Lo,L_Hi)
		pois<-round(rgamma(n, lambda))	
		
		P2<-runif(1,0,1)
		sire2<-rbinom(n,pois,P2)
		meanP2<-mean(sire2)
		
	if(abs(meanP2 - ObsMean)>delta){
		discard[k]<-P2
		k<-k+1
		}else	
	if(abs(meanP2 - ObsMean)<=delta){
		posterior[i]<-P2
		Lambda[l]<-lambda
		i<-i+1
		l<-l+1
		}	
	
}	
	list(posterior = posterior, Lambda = Lambda)
}
