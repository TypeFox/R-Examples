maxnodenum<-function(dendat,h,N,n,d)
{
minim<-matrix(0,d,1)
maxim<-matrix(0,d,1)
i<-1
while (i<=d){
  minim[i]<-min(dendat[,i])  
  maxim[i]<-max(dendat[,i])
  i<-i+1;
}
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)
mindelta<-min(delta)
maxpositive<-ceiling(n*(2*hmax/mindelta)^d)
bigd<-sum(log(N,base=2))
maxnode<-ceiling(bigd*maxpositive)

return(list(maxnode=maxnode,maxpositive=maxpositive));
}
