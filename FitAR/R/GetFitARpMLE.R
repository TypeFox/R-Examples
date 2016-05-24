`GetFitARpMLE` <-
function(z, pvec){
stopifnot(length(z)>0)
if ((length(pvec)==1 && pvec==0) || length(pvec)==0){
    phiHat<-numeric(0)
    Iq<-TRUE
    constantTerm<-mean(z)
    LL<-LoglikelihoodAR(phiHat,z, MeanValue=constantTerm)
    res<-z
    }
else {
	P <- max(pvec)
	stopifnot(length(z)>P)
    ind <- rep(0, P+1)
    ind[pvec] <- NA
    ind[P+1] <- NA
    out<-arima(z, order=c(P, 0, 0), fixed=ind,transform.pars=FALSE)
    res<-resid(out)
    estimates<-coef(out)
    constantTerm<-as.vector(estimates[P+1])
    if (P>0) {
        phiHat<-estimates[1:P]
        Iq<-InvertibleQ(phiHat)
        if (Iq)
            LL<-LoglikelihoodAR(phiHat,z, MeanValue=constantTerm)
        else
            LL<- -1e30
            }
    }
list(loglikelihood=LL, phiHat=phiHat, constantTerm=constantTerm,res=res,pvec=pvec,InvertibleQ=Iq)
}
