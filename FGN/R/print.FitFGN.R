`print.FitFGN` <-
function(x, ...){
H<-x$H
Rsq<-x$Rsq
muHat<-x$muHat
LL<-x$loglikelihood
k<-1
if (!is.null(x$demean)&&x$demean)
    k<-k+1
n<-length(x$res)
aic<- -2*LL+2*k
bic<- -2*LL+log(n)*k
dati<-x$DataTitle
if (!is.null(dati))
    cat(dati,fill=TRUE)
cat(paste("H = ", round(H,3),",  R-sq = ",round(100*Rsq,2),"%", sep=""),fill=TRUE)
cat(paste("length of series =",n, ",  number of parameters =",k),fill=TRUE)
OUTIC<-paste("loglikelihood =",round(LL,3),",  AIC =", round(aic,1),",  BIC = ",round(bic,1))
cat(OUTIC, fill=TRUE)
}

