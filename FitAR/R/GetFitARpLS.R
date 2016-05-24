`GetFitARpLS` <-
function(z, pvec){
stopifnot(length(z)>0)
if (length(pvec)==0 || pvec==0){
    phiHat<-numeric(0)
    constantTerm<-mean(z)
    res<-z-constantTerm
    yX<-matrix(z,ncol=1)
    Iq<-TRUE
    LL<-LoglikelihoodAR(0,z, MeanValue=constantTerm)
    covHat<-numeric(0)
    }
else {
	PMAX<-max(pvec)
	stopifnot(PMAX < length(z))
    Xy <- embed(z, PMAX+1)
    y <- Xy[,1]
    if (length(pvec)==1 && pvec == 1)
        X <- matrix(Xy[,-1], ncol=1)
    else    
        X <- (Xy[,-1])[,pvec]
    yX<-cbind(y,X)
    ans <- lm(y~X)
    res <- resid(ans)
    betaHat <- as.vector(coef(ans)[-1])
    constantTerm <- as.vector(coef(ans)[1])
    phiHat<-numeric(PMAX)
    phiHat[pvec]<-betaHat
    covHat <- vcov(ans)
    }
    Iq<-InvertibleQ(phiHat)
    if (Iq)
        LL<-LoglikelihoodAR(phiHat,z, MeanValue=mean(z))
    else
        LL<- -1e30
list(loglikelihood=LL, phiHat=phiHat, constantTerm=constantTerm, res=res, pvec=pvec,
     InvertibleQ=Iq,  yX=yX, covHat=covHat)
}
