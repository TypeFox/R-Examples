EuCallBS <-
function (Spot,sigma,Strike,r,dr,mT){
	forward<-Spot*exp((r-dr)*mT)
	dplus<-(log(forward/Strike)+0.5*sigma*sigma*mT)/(sigma*sqrt(mT))
	dminus<-(log(forward/Strike)-0.5*sigma*sigma*mT)/(sigma*sqrt(mT))
	mdplus<-dplus*(-1)
	mdminus<-dminus*(-1)
	Ndpl<-pnorm(dplus)
	Ndmin<-pnorm(dminus)
	Nmdpl<-pnorm(mdplus)
	Nmdmin<-pnorm(mdminus)
	V0_BS<-exp(-1*r*mT)*(forward*Ndpl-Strike*Ndmin)
	V0P_BS<-exp(-1*r*mT)*(Strike*Nmdmin-forward*Nmdpl)
	return(V0_BS)
}
