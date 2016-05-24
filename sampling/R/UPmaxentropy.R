"UPmaxentropy" <-function(pik) 
{ 
if(is.data.frame(pik)) 
    if(ncol(pik)>1) stop("pik is not a vector") else pik=unlist(pik)
else if(is.matrix(pik))
       if(ncol(pik)>1) stop("pik is not a vector") else pik=pik[,1]
else  if(is.list(pik))
   if(length(pik)>1) stop("pik is not a vector") else pik=unlist(pik)
n=sum(pik)
n=.as_int(n)
if(n>=2)
{
pik2=pik[pik!=1]
n=sum(pik2)
n=.as_int(n)
piktilde=UPMEpiktildefrompik(pik2) 
w=piktilde/(1-piktilde) 
q=UPMEqfromw(w,n)
s2=UPMEsfromq(q)
s=rep(0,times=length(pik))
s[pik==1]=1
s[pik!=1][s2==1]=1
}
if(n==0) s=rep(0,times=length(pik))
if(n==1) s=as.vector(rmultinom(1, 1,pik)) 
s
}

