m.psamp <-
function(fit){
	B=fit$B1.psamp
	X=unique(fit$matrix.mean)
	n=dim(X)[1]
	p=dim(B)[1]
	nsave=dim(B)[3]
	m.psamp=array(dim=c(n,p,nsave))
	for (i in 1:nsave){
		m.psamp[,,i]=X%*%t(B[,,i])
	}
m.psamp
}
