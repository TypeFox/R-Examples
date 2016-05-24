"FitARMA" <-
function(z,order=c(0,0,0),demean=TRUE,MeanMLEQ=FALSE,pApprox=30,MaxLag=30){
p<-order[1]
d<-order[2]
q<-order[3]
Z<-z
if (d > 0) Z<-diff(z, differences=d)
if (demean) 
    mz<-mean(Z)
else
    mz<-0
y<-Z-mz
pApp<-pApprox
if (q == 0) pApp<-p
ans<-GetFitARMA(y,p,q,pApp)
LL<-ans$loglikelihood
mu<-iter<-0
if (MeanMLEQ && (p>0||q>0)) {
    etol <- MaxIter <- 10
    while(etol> 1e-06 && iter<MaxIter){
        LLPrev<-LL
        iter<-iter+1
        g<-TacvfARMA(ans$phi,ans$theta,pApp)
        coefAR<-PacfDL(g, LinearPredictor=TRUE)$ARCoefficients
        mu<-GetARMeanMLE(y, coefAR)        
        ans<-GetFitARMA(y-mu,p,q,pApp)
        LL<-ans$loglikelihood
        etol<-abs(LL-LLPrev)/LLPrev
        if (ans$convergence != 0) 
            stop("GetARFit returned convergence = ",ans$convergence)           
      }
}
muHat<-mu+mz
phiHat<-ans$phiHat
thetaHat<-ans$thetaHat
if (p>0 || q>0) {
    g<-TacvfARMA(phiHat,thetaHat,pApp)
    coefAR<-PacfDL(g, LinearPredictor=TRUE)$ARCoefficients
    res<-BackcastResidualsAR(y, coefAR, Q=100, demean=FALSE)
    }
else
    res<-Z-muHat
fits<-Z-res
n<-length(res)
sigsq<-sum(res^2)/n
if (p>0 || q>0) 
    covHat<-solve(InformationMatrixARMA(phiHat,thetaHat))/n
else 
    covHat<-numeric(0)
racf<-(acf(res, plot=FALSE, lag.max=MaxLag)$acf)[-1]
LBQ<-LjungBoxTest(res, lag.max=MaxLag, k=order[1]+order[3])
if (d == 0)
    if (q == 0) ModelTitle<-paste("AR(",p,")",sep="")
    if (q > 1) ModelTitle<-paste("ARMA(",p,",",q,")",sep="")
else
    ModelTitle<-paste("ARIMA(",p,",",d,",",q,")",sep="")
out<-list(loglikelihood=ans$loglikelihood,phiHat=phiHat,thetaHat=thetaHat,sigsqHat=sigsq,muHat=muHat,covHat=covHat,
          racf=racf, LjungBoxQ=LBQ,res=res,fits=fits, demean=demean,
          iterationCount=iter,convergence=ans$convergence, MeanMLE=MeanMLEQ, tsp=tsp(z),order=order,call=match.call(),
          DataTitle=attr(z,"title"),ModelTitle=ModelTitle)
class(out)<-"FitARMA"
out
}

