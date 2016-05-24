ABC_P2_gamma <-
function(n,ObsMean, S_Lo, S_Hi, R_Lo, R_Hi, delta,iter){
	posterior<-c()
	discard<-c()
	gamma<-c()
	Shape<-c()
	Rate<-c()
	i<-1
	j<-1
	k<-1
	l<-1
	m<-1

	while(i <= iter){
		dispersion<-runif(1,R_Lo,R_Hi)
		mean<-runif(1,S_Lo,S_Hi)
		gamma<-round(rgamma(n, shape=mean, rate=dispersion))	
		
		P2<-runif(1,0,1)
		sire2<-rbinom(n,gamma,P2)
		meanP2<-mean(sire2)
		
	if(abs(meanP2 - ObsMean)>delta){
		discard[k]<-P2
		k<-k+1
		}else	
	if(abs(meanP2 - ObsMean)<=delta){
		posterior[i]<-P2
		Shape[l]<-mean
		Rate[m]<-dispersion
		i<-i+1
		l<-l+1
		m<-m+1
		}	
	
}	
	list(posterior = posterior, Shape = Shape, Rate = Rate)
}
