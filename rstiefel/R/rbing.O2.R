rbing.O2 <-
function(A,B,a=NULL,E=NULL)
{
  ### assumes B is a diagonal matrix with *decreasing* entries 
  if(is.null(a))
  {
    trA<-A[1,1]+A[2,2]
    lA<-2*sqrt(trA^2/4-A[1,1]*A[2,2]+A[1,2]^2 )
    a<-lA*(B[1,1]-B[2,2])
    E<-diag(2) ; if(A[1,2]!=0)
    {
      E<-cbind(c(.5*(trA+lA)-A[2,2],A[1,2]), c(.5*(trA-lA)-A[2,2],A[1,2]) )
      E[,1]<-E[,1]/sqrt(sum(E[,1]^2)) ; E[,2]<-E[,2]/sqrt(sum(E[,2]^2))
    }
  }

  b<-min(1/a^2,.5) ; beta<- .5-b
  lrmx<- a
  if(beta>0) { lrmx<-lrmx + beta*(log(beta/a)-1) }
  lr<- -Inf
  while(lr< log(runif(1)))
  {
    w<-rbeta(1,.5,b)
    lr<- a*w + beta*log(1-w)   - lrmx
  }
  u<-c(sqrt(w),sqrt(1-w) )*(-1)^rbinom(2,1,.5)
  x1<-E%*%u ; x2<-(x1[2:1]*c(-1,1)*(-1)^rbinom(1,1,.5))
  cbind(x1,x2)
}
