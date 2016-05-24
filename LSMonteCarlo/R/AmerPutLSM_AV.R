AmerPutLSM_AV<-function (Spot=1, sigma=0.2, n=1000, m=365, Strike=1.1, r=0.06, dr=0.0, mT=1) {
	Z<-matrix(rnorm(m*n, mean= 0, sd=1), nrow=n, ncol=m)
	aZ<-Z*(-1)
	ZZ<-rbind(Z, aZ)
	GBM<-matrix(NA, nrow=n*2, ncol=m)	
	for(i in 1:(n*2)) {
	GBM[i,]<-Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*ZZ[i,])))
	}
	X<-ifelse(GBM<Strike,GBM,NA)
	CFL<-matrix(pmax(0,Strike-GBM), nrow=n*2, ncol=m)
	Xsh<-X[,-m]
	X2sh<-Xsh*Xsh
	Y1<-CFL*exp(-1*r*(mT/m))
	Y2<-cbind((matrix(NA, nrow=n*2, ncol=m-1)), Y1[,m])
	CV<-matrix(NA, nrow=n*2, ncol=m-1)
	try(for(i in (m-1):1) {
	reg1<-lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
	CV[,i]<-(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
	CV[,i]<-(ifelse(is.na(CV[,i]),0,CV[,i]))
	Y2[,i]<-ifelse(CFL[,i]>CV[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(mT/m)))
	}
	, silent = TRUE)
	CV<-ifelse(is.na(CV),0,CV)
	CVp<-cbind(CV, (matrix(0, nrow=n*2, ncol=1)))
	POF<-ifelse(CVp>CFL,0,CFL)
	FPOF<-firstValueRow(POF)
	dFPOF<-matrix(NA, nrow=n*2, ncol=m)
	for(i in 1:m) {
		dFPOF[,i]<-FPOF[,i]*exp(-1*mT/m*r*i)
		}
	PRICE<-mean(rowSums(dFPOF))
	res<- list(price=(PRICE), Spot, Strike, sigma, n, n, m, r, dr, mT)
	class(res)<-"AmerPutAV"
	return(res)
}
