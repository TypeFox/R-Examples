"pcopula1" <-
function(theta,delta,psi,phi,ivpsi,ivphi,s,t)
{n<-min(length(s),length(t));
	u<-ivphi(theta,delta,exp(-ivpsi(delta,s)[(n+2):(2*n+1)]))[(2+n+1):(2*n+2)];
	v<-ivphi(theta,delta,exp(-ivpsi(delta,t)[(n+2):(2*n+1)]))[(2+n+1):(2*n+2)];
	suma<-u+v;
	pc<-psi(delta,-log(phi(theta,delta,suma)[(2+n+1):(2*n+2)]))[(n+2):(2*n+1)]
	resu<-pc
}

