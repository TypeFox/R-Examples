`print.FitAR` <-
function(x, ...){
LL<-x$loglikelihood
k<-length(x$pvec)
if (!is.null(x$demean)&&x$demean)
    k<-k+1
n<-length(x$res)
aic<- -2*LL+2*k
bic<- -2*LL+log(n)*k
subQ<-x$SubsetQ
if (subQ) {
    lags<-x$pvec
    P<-max(lags)
    ubic<-bic + 2*lchoose(P, k)
    }
dati<-x$DataTitle
if (!is.null(dati))
    cat(dati,fill=TRUE)
modti<-x$ModelTitle
    if (x$FitMethod=="LS")
        modti<-paste("AR(",max(x$pvec),"). LS Fit.",sep="")
    else
        modti<-paste("AR(",max(x$pvec),"). MLE.",sep="")
if (x$ARModel=="ARz") 
    if (x$MeanMLE)
        modti<-paste(modti, " Mean estimated using MLE")
    else
        modti<-paste(modti, " Mean estimated using the sample mean")
cat(modti,fill=TRUE)
cat(paste("length of series =",n, ",  number of parameters =",k),fill=TRUE)
OUTIC<-paste("loglikelihood =",round(LL,3),",  AIC =", round(aic,1),",  BIC = ",round(bic,1))
if (subQ) 
    OUTIC<-paste(OUTIC, ", UBIC = ", round(ubic,1))
cat(OUTIC, fill=TRUE)
}

