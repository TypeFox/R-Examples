Multiple.Testing <-
function(s1,s2,m,p,D,delta,BCS,pho,K,alpha,beta){

L=ceiling(m/BCS);length=m;
delta.vec=rep(0,length);
delta.vec[1:D]=delta;

epsilon<-function(length,L,m,pho){
u.vec=rnorm(length)
b=rnorm(L)
b.vec=0
for (i in 1:L)
{
b.vec[((i-1)*100+1):(i*100)]=b[i]
}
b.vec=b.vec[1:m]
epsilon.vec=u.vec*sqrt(1-pho)+b.vec*sqrt(pho);
return(epsilon.vec)
}
# generate K replicates error vector
epsilon.maximun=rep(0,K)
for (k in 1:K)
{
epsilon.maximun[k]=max(epsilon(length,L,m,pho))
}

c_alpha=as.numeric(quantile(epsilon.maximun, min(((1-alpha)*K+1)/K,1)))

h<-function(n)
{
h=rep(0,K)
h.max=rep(0,K)
for(k in 1:K)
{
h.max[k]=max(epsilon(length,L,m,pho)+delta.vec*sqrt(n*p*(1-p)))
h[k]=ifelse(h.max[k]>c_alpha,1,0)
}
hh=sum(h)/K -1 + beta
return(hh)
}

if (h(s1) < 0 & h(s2)>0)
{
while(h(s1) < 0 & h(s2)>0 & (abs(s1-s2)>1 | abs(h(s2))>1)) 
{ 
s3=(s1+s2)/2; 
if(h(s3)>0) s2=s3 else s1=s3 
}
}
return(list("s1"=s1,"s2"=s2,"h(s1)"=h(s1),"h(s2)"=h(s2)))
}
