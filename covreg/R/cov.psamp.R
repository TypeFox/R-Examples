cov.psamp <-
function(fit){
	A.psamp=fit$A.psamp
	B.psamp=fit$B2.psamp
	nsave=dim(A.psamp)[3]
	p=dim(A.psamp)[1]
	R=dim(B.psamp)[3]
	X=unique(fit$matrix.cov)
	n=dim(X)[1]
	s.psamp=array(dim=c(n,p,p,nsave))
	for (iter in 1:nsave){
		for (i in 1:n){
			ss=A.psamp[,,iter]
			for (r in 1:R){
				ss=ss+B.psamp[,,r,iter]%*%X[i,]%*%t(X[i,])%*%t(B.psamp[,,r,iter])
			}
			s.psamp[i,,,iter]=ss
		}
	}
s.psamp
}
