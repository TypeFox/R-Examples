`gBox` <-
function(model,lags=1:20,x=eval(model$call$x), method=c('squared','absolute')[1], plot=TRUE){

# Carry out the generalized Box test for residual autocorrelation in the residuals
# from a fitted GARCH model.
#

gBox.test <-
function(model,lag,x)
{
p=model$order[1]
q=model$order[2]
start=1+q
beta=model$coef[(start+1):(start+p)]
r=x
epsilon=residuals(model)
M=NULL
i=1
while(i<=q) {M=cbind(M,zlag(r^2,i));i=i+1}
i=1
h=model$fit[,1]
h=h*h
while(i<=p) {M=cbind(M,zlag(h,i));i=i+1}
sigma2=mean(r^2)
M[is.na(M)]=sigma2
if(p>0) M=filter(M,filter=beta,method='recursive',sides=1,init=rep(sigma2,
length(beta)))
if(p>0) M=cbind(1/(1-sum(beta)),M) else M=cbind(1,M)
M=apply(M,2,function(x){x/h})
M[is.na(M)]=0
E=NULL
for (i in 1:lag) {E=cbind(E,zlag(epsilon^2,i)-mean(epsilon^2,na.rm=TRUE))}
E[is.na(E)]=0
J=t(E)%*%M/length(h)
k=kurtosis(epsilon,na.rm=TRUE)+2
Lambda=2*solve(t(M)%*%M)
H=-J%*%Lambda%*%t(J)/(2*k)*length(h)
diag(H)=1+diag(H)
acf1=matrix(as.vector(acf(epsilon^2,lag.max=lag,na.action=na.omit,plot=FALSE)$acf),ncol=1)
stat=as.numeric(t(acf1)%*%solve(H)%*%acf1*length(h))
p.value=1-pchisq(stat,lag)
return(invisible(list(p.value=p.value,H=H,Lambda=Lambda,J=J,
stat=stat,acf=acf1,n=length(h))))
}
gBox1.test<-
function(model,lag,x)
{
p=model$order[1]
q=model$order[2]
start=1+q
beta=model$coef[(start+1):(start+p)]
r=x
epsilon=residuals(model)
M=NULL
i=1
while(i<=q) {M=cbind(M,zlag(r^2,i));i=i+1}
i=1
h=model$fit[,1]
h=h*h
while(i<=p) {M=cbind(M,zlag(h,i));i=i+1}
sigma2=mean(r^2)
M[is.na(M)]=sigma2
if(p>0) M=filter(M,filter=beta,method='recursive',sides=1,init=rep(sigma2,
length(beta)))
if(p>0) M=cbind(1/(1-sum(beta)),M) else M=cbind(1,M)
M=apply(M,2,function(x){x/h})
M[is.na(M)]=0
E=NULL
tau=mean(abs(epsilon),na.rm=TRUE)
nu=mean(abs(epsilon)^3,na.rm=TRUE)
for (i in 1:lag) {E=cbind(E,zlag(abs(epsilon),i)-tau)}
E[is.na(E)]=0
J=t(E)%*%M/length(h)
k=kurtosis(epsilon,na.rm=TRUE)+2
Lambda=2*solve(t(M)%*%M)*k/2
H=J%*%Lambda%*%t(J)*length(h)*(k/2*tau^2/4-tau*(nu-tau)/2)/(1-tau^2)^2
diag(H)=1+diag(H)
acf1=matrix(as.vector(acf(abs(epsilon),lag.max=lag,na.action=na.omit,plot=FALSE)$acf),ncol=1)
stat=as.numeric(t(acf1)%*%solve(H)%*%acf1*length(h))
p.value=1-pchisq(stat,lag)
return(invisible(list(p.value=p.value,H=H,Lambda=Lambda,J=J,
stat=stat,acf=acf1,n=length(h))))
}

pv=NULL
for (lag in lags) {
new=switch(method, squared=gBox.test(model,lag,x)$p.value,
absolute=gBox1.test(model,lag,x)$p.value)
pv=c(pv,new)
}
if (plot) {plot(x=1:20,y=pv,ylim=c(0,1),xlab='lag',ylab='p-value',main='')
abline(h=0.05,lty=2)}
invisible(list(pvalue=pv,lags=lags,method=method,x=x))
}

