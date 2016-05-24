`changedimension_covariates` <-
function(X,d,r,n){
	covariates2 = array(NA,dim=c(d,r,n))
	seq1=seq(1,(d*n)+(-n+1),by=n)
	seq2=seq(n,(d*n),by=n)

	for(i in 1:d){
		covariates2[i,,] = t(X[seq1[i]:seq2[i],])}
	return(covariates2)	
}

