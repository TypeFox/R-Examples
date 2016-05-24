"UPminimalsupport" <-
function(pik)
{
if(any(is.na(pik))) stop("there are missing values in the pik vector")
basicsplit<-function(pik)
{
N=length(pik)
n=sum(pik)
A=(1:N)[pik==0]
B=(1:N)[pik==1]
C=setdiff(setdiff(1:N,A),B)
D=C[sample(length(C), round(n-length(B)))]
s1v=rep(0,times=N)
s1v[c(B,D)]=1
alpha=min(1-max(pik[setdiff(C,D)]),min(pik[D]))
pikb= (pik-alpha*s1v)/(1-alpha)
if(runif(1,0,1)<alpha) s=s1v else s=pikb
s
}
is.a.sample<-function(s,EPS=sqrt(.Machine$double.eps)) if(sum(abs(s-round(s)))<EPS) TRUE else FALSE
while(!is.a.sample(pik))pik=basicsplit(pik) 
round(pik)
}

