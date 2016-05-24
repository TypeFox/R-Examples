rmf.vector <-
function(kmu)
{
  #simulate from the vector mf distribution as described in Wood(1994)
  kap<-sqrt(sum(kmu^2)) ; mu<-kmu/kap ; m<-length(mu)
  if(kap==0){ u<-rnorm(length(kmu)) ; u<-matrix(u/sqrt(sum(u^2)),m,1) }
  if(kap>0)
  {
    if(m==1){ u<- (-1)^rbinom( 1,1,1/(1+exp(2*kap*mu))) }
    if(m>1)
    {
      W<-rW(kap,m)
      V<-rnorm(m-1) ;  V<-V/sqrt(sum(V^2))
      x<-c((1-W^2)^.5*t(V),W)
      u<-cbind( NullC(mu),mu)%*%x
    }
  }
  u
}
