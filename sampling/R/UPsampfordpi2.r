UPsampfordpi2<-function(pik) 
{
n=sum(pik)
n=.as_int(n)
if(n<2) stop("the sample size<2")
N=length(pik)
p=pik/n
pikl=matrix(0,N,N)
Lm=rep(0, n)
lambda=p/(1-n*p)
Lm[1]=1
if(n>=2)
for (i in 2:n) {
for (r in 1:(i-1)) Lm[i]=Lm[i]+((-1)^(r-1))*sum(lambda^r)*Lm[i-r]
Lm[i]=Lm[i]/(i - 1)
}
if(any(Lm<0)) stop("it is not possible to compute pik2 for this example")
t1=(n + 1) - (1:n)
Kn=1/sum(t1*Lm/n^t1)
Lm2=rep(0, n - 1)
t2=(1:(n - 1))
t3=n - t2 
for (i in 2:N) {
for (j in 1:(i - 1)) {
Lm2[1]=1
Lm2[2]=Lm[2] - (lambda[i] + lambda[j])
if(n>3)
for (m in 3:(n - 1)) {
Lm2[m]=Lm[m] - (lambda[i] + lambda[j]) * Lm2[m -1] - lambda[i] * lambda[j] * Lm2[m - 2]
}
pikl[i, j]=Kn * lambda[i] * lambda[j] * sum((t2+1-n*(p[i] + p[j]))*Lm2[t3]/n^(t2 - 1))
pikl[j, i]=pikl[i, j]
}
pikl[i, i]=pik[i] 
}
pikl[1, 1]=pik[1] 
pikl
}
