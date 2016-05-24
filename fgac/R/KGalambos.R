"KGalambos" <-
function(u,v,delta)
{
	utilde<--log(u);
	vtilde<--log(v);
	KGalambos<-u*v*exp((utilde^(-delta)+vtilde^(-delta))^(-1/delta));
	resu<-matrix(c(u,v,delta,KGalambos),nrow=1)
	
}

