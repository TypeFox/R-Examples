"UPtille" <-
function(pik,eps=1e-6)
{
if(any(is.na(pik))) stop("there are missing values in the pik vector")
n=sum(pik)
n=.as_int(n)
list = pik > eps & pik < 1 - eps
pikb = pik[list]
N = length(pikb)
s=pik
if(N<1) stop("the pik vector has all elements outside of the range [eps,1-eps]")
else 
{
n=sum(pikb)
sb=rep(1,N)
b=rep(1,N)
for(i in 1:(N-n))
	{a=inclusionprobabilities(pikb,N-i)
	v=1-a/b
	b=a
	p=v*sb
	p=cumsum(p)
	u=runif(1)
        for(j in 1:length(p))
		if(u<p[j]) break
        sb[j] = 0
	}
s[list]=sb
}
s
}


