"UPtillepi2" <-
function(pik,eps=1e-6)
{

if(any(is.na(pik))) warning("there are missing values in the pik vector")
n=sum(pik)
n=.as_int(n)
list = pik > eps & pik < 1 - eps
pikb = pik[list]
N = length(pikb)
#ppf=pik%*%t(pik)
ppf=matrix(0,length(pik),length(pik))
if(N<1) stop("the pik vector has all elements outside of the range [eps,1-eps]")
else 
{
n=sum(pikb)
if(N>n)
{
UN=rep(1,N)
b=rep(1,N)
pp=1
for(i in 1:(N-n))
	{
	a=inclusionprobabilities(pikb,N-i)
	vv=1-a/b
	b=a
	d=vv %*% t(UN)
	pp=pp*(1-d-t(d))
	}
diag(pp)=pikb
ppf[list,list]=pp
}
}
ppf
}

