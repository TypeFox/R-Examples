BICqSelect <-
function(logL, n, level=0.99, mSize=1:length(logL), mComplex=function(k) k)
{
#AICa = -2logL(k) + a*C(k)
#BICq = -2logL(k) + C(k)*[log(n)-2log{q/(1-q)}]
#C(k) is the model complexity or degree of freedom, usually C(k)=k.
#logL: log-likelihood
#n: sample size
#level: confidence level for controlling overfitting
#multiple level can be used, such as level=c(0.95, 0.99).
#mSize: the set of model sizes {k1, k2, ..., k_P}, usually mSize=1:P.
#Model complex: 
mC <- mComplex(mSize) 
#The number of candidate models:
P<- length(logL)
stopifnot(P>2)
stopifnot(length(unique(mSize))==P) #Model sizes are not unique.
if (!all(mSize[2:P]-mSize[1:(P-1)]>=0) && P>1) stop("Model sizes are not ordered.")
#The tuning parameter a for AICa:
aChosen<- qchisq(level,1) 
## Ranges of a for AICa 
a12<- matrix(rep(NA,P*2),ncol=2)
a12[1,]<- 2*c(max((logL[2:P]-logL[1])/(mC[2:P]-mC[1])),Inf)
a12[P,]<- 2*c(0, min((logL[1:(P-1)]-logL[P])/(mC[1:(P-1)]-mC[P])))
for (k in 2:(P-1)){
    i1<-1:(k-1) 
    i2<-(k+1):P
    a12[k,]<- 2*c(max((logL[i2]-logL[k])/(mC[i2]-mC[k])), 
               min((logL[i1]-logL[k])/(mC[i1]-mC[k])))
    }
## Ranges of q for BICq
q12 <- 1/(1+exp(a12/2)/sqrt(n))
q12 <- q12[,2:1]
## Select the best model
ks<- rep(NA, length(aChosen))
as<- matrix(rep(NA, 2*length(aChosen)),ncol=2)
qs<- matrix(rep(NA, 2*length(aChosen)),ncol=2)
for (i in 1:length(aChosen)){
    indx<- which(a12[,1]<=aChosen[i] & a12[,2]>=aChosen[i])
    ks[i]<- mSize[indx]
    as[i,]<- a12[indx,]
    qs[i,]<- q12[indx,]
    }
kaqChosen<- data.frame(k=ks, a=as, q=qs, level=level)
aqTable<- data.frame(k=mSize, a=a12, q=q12, ms=a12[,1]<=a12[,2])
colnames(aqTable)[6]="a1<=a2"
list(kHat=kaqChosen, table=aqTable)
}

