`FitFGN` <-
function(z,demean=TRUE,MeanMLEQ=FALSE,lag.max="default"){
n<-length(z)
stopifnot(n>0)
MaxIter<-10 #max iterations for exact MLE
ztsp<-tsp(z)
if (lag.max=="default")
    MaxLag=ceiling(min(length(z)/4,min(max(length(z)/4,30),100)))
else
    MaxLag=lag.max
MaxIter<-10
indMeanQ <- demean || MeanMLEQ
if (indMeanQ)
    mz <- mean(z) 
else
    mz <- 0
y <- z - mz
ans<-GetFitFGN(y,MeanZeroQ=TRUE)
LL<-ans$loglikelihood
etol <- 1
mu<- iter <- 0
if (MeanMLEQ)
    while(etol> 1e-06 && iter<MaxIter){
        LLPrev<-LL
        iter<-iter+1
        r<-acvfFGN(ans$H, n-1)
        mu<-TrenchMean(r, y)
        ans<-GetFitFGN(y-mu,MeanZeroQ=TRUE)
        LL<-ans$loglikelihood
        etol<-abs(LL-LLPrev)/LLPrev
        if (ans$convergence != 0) 
            stop("GetARFit returned convergence = ",ans$convergence)           
      }
muHat<-mu+mz
H<-ans$H
SEH<-sqrt((c(0.13,0.24,0.31,0.37,0.40,0.43,0.45,0.47,0.47, 0.47)[max(1,round(10*H))])/n)
rH<-var(z)*acvfFGN(H, n-1)
SEmu<-sqrt((rH[1]+2*sum(rH[-1]*(n-(1:(n-1))))/n)/n)
res<-DLResiduals(rH, y)
racf<-(acf(res, plot=FALSE, lag.max=MaxLag)$acf)[-1]
L0<-10
lags<-L0:MaxLag
QQ<-n*(n+2)*cumsum((racf^2)/(n-(1:MaxLag)))[lags]
pv<-1-pchisq(QQ,lags-1)
QQ<-round(QQ,2)
a<-matrix(c(lags,QQ,pv),ncol=3)
dimnames(a)<-list(rep("",length(QQ)),c("m","Qm", "pvalue"))
LBQ<-a
if (H < 0.75)
    rH2<-acvfFGN(H, 10^3)
else
    rH2<-acvfFGN(H, 10^4)
Rsq<-1-DLAcfToAR(rH2[-1])[-1+length(rH2),3]
sigsq<-var(z)*(1-Rsq)
ans<-list(loglikelihood=ans$loglikelihood, H=H, SEH=SEH, sigsqHat=sigsq,muHat=muHat,SEmu=SEmu, Rsq=Rsq,
          LjungBoxQ=LBQ,res=res,demean=demean, 
          iterationCount=iter, convergence=ans$convergence, MeanMLEQ=MeanMLEQ, z=z, tsp=ztsp, call=match.call(),
          DataTitle=attr(z,"title"))
class(ans)<-"FitFGN"
ans
}
