`summary.FitFGN` <-
function(object, ...){
H<-object$H
Rsq<-object$Rsq
muHat<-object$muHat
sigsq<-object$sigsq
LL<-object$loglikelihood
k<-1
if (!is.null(object$demean)&&object$demean)
    k<-k+1
n<-length(object$res)
aic<- -2*LL+2*k
bic<- -2*LL+log(n)*k
dati<-object$DataTitle
if (!is.null(dati))
    cat(dati,fill=TRUE)
cat(paste("H = ", round(H,3),",  R-sq = ",round(100*Rsq,2),"%", sep=""),fill=TRUE)
cat(paste("length of series =",n, ",  number of parameters =",k),fill=TRUE)
cat(paste("mean = ", muHat, "RMSE = ",sqrt(sigsq)),fill=TRUE)
OUTIC<-paste("loglikelihood =",round(LL,3),",  AIC =", round(aic,1),",  BIC = ",round(bic,1))
cat(OUTIC, fill=TRUE)
LB<-object$LjungBoxQ
cat(LB[nrow(LB),], fill=TRUE)
}

