dw.meanvar<-function(q,beta, M=1000)
{
idx<-1:M
m.res<-sum(as.numeric(lapply(idx,function(x,q,beta) q^(x^(beta)),q=q,beta=beta)))
x2<-sum(as.numeric(lapply(idx,function(x,q,beta) x*q^(x^(beta)),q=q,beta=beta)))
v.res<-2*x2-m.res-m.res^2
return(list(mean=m.res,var=v.res))
}