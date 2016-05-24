`summary.FitAR` <-
function(object, ...){
RSQ<- 1- (object$sigsqHat)/var(object$z)
LL<-object$loglikelihood
k<-length(object$pvec)
if (!is.null(object$demean)&&object$demean)
    k<-k+1
n<-length(object$res)
aic<- -2*LL+2*k
bic<- -2*LL+log(n)*k
subQ<-object$SubsetQ
if (subQ) {
    lags<-object$pvec
    P<-max(lags)
    ubic<-bic + 2*lchoose(P, k)
    }
dati<-object$DataTitle
if (!is.null(dati))
    cat(dati,fill=TRUE)
modti<-object$ModelTitle
    if (object$FitMethod=="LS")
        modti<-paste("AR(",max(object$pvec),"). LS Fit.",sep="")
    else
        modti<-paste("AR(",max(object$pvec),"). MLE.",sep="")
if (object$FitMethod=="MLE" && object$MeanMLE) 
    modti<-paste(modti, " With mean MLE.")    
cat(modti,fill=TRUE)
cat(paste("length of series =",n, ",  number of parameters =",k),fill=TRUE)
OUTIC<-paste("loglikelihood =",round(LL,3),",  aic =", round(aic,1),",  bic = ",round(bic,1))
if (subQ) 
    OUTIC<-paste(OUTIC, ", UBIC = ", round(ubic,1))
cat(OUTIC, fill=TRUE)
cat(paste("series mean = ", object$muHat, ",  rmse = ", sqrt(object$sigsqHat),
   ",  R^2 = ", 100*round(RSQ,5), "%", sep=""), fill=TRUE)
invisible()
}

