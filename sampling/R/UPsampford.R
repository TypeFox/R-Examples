UPsampford<-function(pik,eps=1e-6,max_iter=500)
{
if(any(is.na(pik))) stop("there are missing values in the pik vector")
n=sum(pik)
n=.as_int(n)
list= pik>eps & pik < 1-eps
pikb=pik[list]
n=sum(pikb)
N=length(pikb)
s=pik
if(N<1) stop("the pik vector has all elements outside of the range [eps,1-eps]")
else 
{
sb=rep(2,N)
y=pikb/(1-pikb)/sum(pikb/(1-pikb))
step=0
while(sum(sb<=1)!=N & step<=max_iter)
      {
	sb=as.vector(rmultinom(1,1,pikb/sum(pikb))+rmultinom(1,.as_int(n-1),y))
      step=step+1
       }
if(sum(sb<=1)==N) s[list]=sb
else stop("Too many iterations. The algorithm was stopped.")
}
s
}
