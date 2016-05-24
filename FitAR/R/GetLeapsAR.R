GetLeapsAR <-
function(z, lag.max=15, Criterion="UBIC", Best=3, Candidates=5, t="default", ExactQ=FALSE){
stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0)
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
stopifnot(is.wholenumber(lag.max))
method<-Criterion
if (is.na(pmatch(method,c("UBIC","AIC","BIC","EBIC","BICq","GIC"))))
    method<-"UBIC"
if (ExactQ){
    stop("Sorry this option is not available yet!")
    M<-2^lag.max
    logL <- numeric(M)
    logL[1] <- GetFitARpMLE(z, pvec=0)$loglikelihood
    for (i in 1:(M-1)){
        ind<-as.logical(rev(toBinary(i, lag.max)))
        pvec <- (1:lag.max)[ind]
        out<-GetFitARpMLE(z, pvec=pvec)
        logL[i+1] <- out$loglikelihood
    }
}
#set tuning parameter
P<-0.01
Q<-0.25
G<-1
if (method=="EBIC"  && t!="default")  G <- t
if (method=="QBIC"  && t!="default")  Q <- t
if (method=="GIC"   && t!="default")  P <- t
if (P>=0.25 || P<=0)
    stop("error: GIC tuning parameter invalid")
#level <- 1 - P
if (Q<=0 || Q>=1)
    stop("error: BICq tuning parameter invalid")
#GIC/BICq treated as a special case
pvec <- 1:lag.max
LagRange<-1:lag.max
n <- length(z)-lag.max
ind <- (lag.max+1):length(z)
y<-z[ind]
X<-matrix(rep(0,n*lag.max), nrow=n, ncol=lag.max)
for (i in 1:lag.max)
    X[,i] <- z[ind-pvec[i]]
outLeaps <- leaps(y=y,x=X,nbest=1,method="r2", strictly.compatible=FALSE)
k <- outLeaps$size
#
#this defines an approximate likelihood approach
TotSS <- sum((y-mean(y))^2)
RSS <- TotSS*(1-outLeaps$r2)
LogL <- (-n/2)*log(RSS/n)
if (method=="AIC")
    ic<- -2*LogL + 2*k
if (method=="BIC")
    ic<- -2*LogL + log(n)*k
if (method=="UBIC")
    ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
if (method=="EBIC")
    ic<- -2*LogL + log(n)*(1+LagRange)+2*G*lchoose(lag.max, k)
if (method=="BICq")
    ic<- -2*LogL + log(n)*(1+LagRange)-2*(LagRange*log(Q)+(lag.max+1-k)*log(1-Q))
if (method=="GIC")
    ic<- -2*LogL + k*qchisq(p=(1+sqrt(1-4*P))/2, df=1) 
indBest<-order(ic)
#extra step needed because leaps does not include null model
#Very important: refit with exact MLE
LogL<-numeric(Candidates+1)
#LogL[1]<-GetFitARpLS(z, 0)$loglikelihood #null model included here
LogL <-GetFitAR(z, 0)$loglikelihood #null model included here
for (i in 1:Candidates)
       LogL[i+1]<-GetFitAR(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
       #LogL[i+1]<-GetFitARpLS(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
k<-c(1,k[indBest[1:Candidates]])
if (method=="AIC")
    ic<- -2*LogL + 2*k
if (method=="BIC")
    ic<- -2*LogL + log(n)*k
if (method=="EBIC")
    ic<- -2*LogL + log(n)*k + 2*G*lchoose(lag.max+1, k)
if (method=="UBIC")
    ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
if (method=="BICq")
    ic<- -2*LogL + k*(log(n) - 2*log(Q/(1-Q)))
if (method=="GIC")
    ic<- -2*LogL + k*qchisq(p=(1+sqrt(1-4*P))/2, df=1)
icBest<-order(ic)[1:Best] #best models exact
m<-as.list(numeric(Best))
for (i in 1:Best){
    ind<-icBest[i]
    if (indBest[ind] ==1)
        p<-0
    else
        p<-pvec[outLeaps$which[indBest[ind-1],]]
    if (method == "AIC")
        m[[i]] <-list(p=p, AIC=ic[ind])
    if (method == "BIC")
        m[[i]] <-list(p=p, BIC=ic[ind])
    if (method == "UBIC")
        m[[i]] <-list(p=p, UBIC=ic[ind])
    if (method == "EBIC")
        m[[i]] <-list(p=p, EBIC=ic[ind], g=G)
    if (method == "BICq")
        m[[i]] <-list(p=p, BICq=ic[ind], q=Q)
    if (method == "GIC")
        m[[i]] <-list(p=p, GIC=ic[ind], p=P)
    }
class(m)<-"Selectmodel"
attr(m,"model")<-"ARp"
if (Best>1)
    ans<-m
else {
    ans <- m[[1]]$p
    if (length(ans)==0)
        ans<-0
    }
ans
}
