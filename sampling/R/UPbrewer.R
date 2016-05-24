"UPbrewer" <- function(pik, eps = 1e-06)
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
sb=rep(0,N)
n=sum(pikb)
for (i in 1:n) {
        a = sum(pikb*sb)
        p = (1-sb)*pikb*((n-a)-pikb)/((n-a)-pikb*(n-i+1))
        p = p/sum(p)
        p = cumsum(p)
        u=runif(1)
        for(j in 1:length(p))
		if(u<p[j]) break
        sb[j] = 1
}
s[list]=sb
}
s
}



