QuantoAmerPutLSM_AV <-
function(Spot=1, sigma=0.2, n=1000, m=365, Strike=1.1, r=0.06, dr=0.0, mT=1, Spot2=1, sigma2=0.2, r2=0.0, dr2=0.0, rho=0){
	covmat <- matrix(c(1,rho,rho,1), ncol=2)
	GBM1a<-matrix(NA, nrow=n, ncol=m)
	GBM1b<-matrix(NA, nrow=n, ncol=m)
	GBM2a<-matrix(NA, nrow=n, ncol=m)
	GBM2b<-matrix(NA, nrow=n, ncol=m)
	GBMD<-matrix(NA, nrow=m, ncol=2)
	for(i in 1:n) {
	GBMDa<-rmvnorm(m, mean=c(0,0), sigma=covmat)
	GBMDb<-GBMDa*(-1)
	GBM1a[i,]<-Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*GBMDa[,1])))
	GBM1b[i,]<-Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*GBMDb[,1])))
	GBM2a[i,]<-Spot2*exp(cumsum(((r2-dr2)*(mT/m)-0.5*sigma2*sigma2*(mT/m))+(sigma2*(sqrt(mT/m))*GBMDa[,2])))
	GBM2b[i,]<-Spot2*exp(cumsum(((r2-dr2)*(mT/m)-0.5*sigma2*sigma2*(mT/m))+(sigma2*(sqrt(mT/m))*GBMDb[,2])))
	}
	GBM1<-rbind(GBM1a,GBM1b)
	GBM2<-rbind(GBM2a,GBM2b)
	X<-ifelse(GBM1<Strike,GBM1,NA)
	CFL1<-matrix(pmax(0,Strike-GBM1), nrow=n*2, ncol=m)
	CFL<-CFL1*GBM2
	Xsh<-X[,-m]
	X2sh<-Xsh*Xsh
	G<-ifelse(GBM1<Strike,GBM2,NA)
	Gsh<-G[,-m]
	G2sh<-Gsh*Gsh
	H<-Xsh*Gsh
	Y1<-CFL*exp(-1*r*(mT/m))
	Y2<-cbind((matrix(NA, nrow=n*2, ncol=m-1)), Y1[,m])
	CV<-matrix(NA, nrow=n*2, ncol=m-1)
	try(for(i in (m-1):1) {
	reg1<-lm(Y2[,i+1]~Xsh[,i]+X2sh[,i]+Gsh[,i]+G2sh[,i]+H[,i])
	mat1<-cbind(matrix(reg1$coefficients)[1,1], matrix(reg1$coefficients)[2,1], matrix(reg1$coefficients)[3,1], matrix(reg1$coefficients)[4,1], matrix(reg1$coefficients)[5,1], matrix(reg1$coefficients)[6,1])
	mat2<-(ifelse(is.na(mat1),0,mat1))
	CV[,i]<-(mat2[1,1])+((mat2[1,2])*Xsh[,i])+((mat2[1,3])*X2sh[,i])+((mat2[1,4])*Gsh[,i])+((mat2[1,5])*G2sh[,i])+((mat2[1,6])*H[,i])
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
	res<- list(price=(PRICE), Spot, Strike, sigma, n, n, m, r, dr, mT, Spot2, sigma2, r2, dr2, rho)
	class(res)<-"QuantoAmerPut_AV"
	return(res)
}
