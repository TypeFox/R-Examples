AmerPutLSMPriceSurf <-
function(Spot=1, vols=(seq(0.1,2,0.1)), n=1000, m=365, strikes=(seq(0.5,2.5,0.1)), r=0.06, dr=0.0, mT=1) {
	volm<-matrix(vols, ncol=1, nrow=(length(vols)))
	strikesm<-matrix(strikes, ncol=1, nrow=(length(strikes)))
	prss1<-matrix(NA, ncol=(length(strikes)), nrow=(length(vols)))
	for(i in 1:(length(vols))) {
		for(j in 1:(length(strikes))) {
			prss1[i,j]<-AmerPutLSM(Spot,volm[i,1],n,m,strikesm[j,1],r,dr,mT)$price
			}
		}
	colnames(prss1)<-strikesm
	rownames(prss1)<-volm
	class(prss1)<-"PriceSurface"
	return(prss1)
}
