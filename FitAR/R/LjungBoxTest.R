`LjungBoxTest` <-
function(res, k=0, lag.max=30, StartLag=1, SquaredQ=FALSE){
stopifnot(k>=0, StartLag>=1, lag.max>=StartLag)
n<-length(res)
L0<-StartLag
if (SquaredQ) {
    z<-(res-mean(res))^2
    kpar<-0
    }
else {
    z<-res
    kpar<-k
}
ra<-(acf(z, lag.max=lag.max, plot=FALSE)$acf)[-1]
lags<-L0:lag.max
QQ<-n*(n+2)*cumsum((ra^2)/(n-(1:lag.max)))[lags]
df <- ifelse(lags-kpar > 0, lags-kpar, 1)
pv<-1-pchisq(QQ,df)
QQ<-round(QQ,2)
a<-matrix(c(lags,QQ,pv),ncol=3)
dimnames(a)<-list(rep("",length(QQ)),c("m","Qm", "pvalue"))
a
}

