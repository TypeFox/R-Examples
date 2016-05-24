"pcopula2" <-
function(theta,delta,psi,v1,ivpsi,v2,s,t)
{if(missing(v1) | missing(v2)){v1<-1;v2<-1};
n<-min(length(s),length(t));
pc<-psi(theta,-log(KGalambos(exp(-ivpsi(theta,s)[(n+2):(2*n+1)]),exp(-ivpsi(theta,t)[(n+2):(2*n+1)]),delta)[(2+2*n):(3*n+1)]))[(n+2):(2*n+1)]
resu<-pc
}

