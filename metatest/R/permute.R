permute <-
function(b,y,x,v,residual,xx,tval,npermut,itar,treshold,maxiter){

	n=length(y)
	p=dim(x)[2]
	tp=c()
		
	for (i in 1:npermut) {
		
		a=sample((1:n))
		bootsamp=rbinom(n,1,0.5)
		for (i in 1:(n)) {if(bootsamp[i]==0.0) bootsamp[i]=-1}
		
		if(p>0) yp=x %*% b +  residual[a]*bootsamp             # multiply by +/- 1 in order to get adequate tests for intercept only model
		else yp=residual[a]*bootsamp             # multiply by +/- 1 in order to get adequate tests for intercept only model
		
		sigmasq=estim(yp,xx,v[a],treshold,maxiter)[1]  # v's are also permuted
		w=diag(1/(v[a] + sigmasq))
		covmat= solve(t(xx) %*% w %*% xx)
		bp=(covmat %*% t(xx) %*% w %*% yp)[itar]
		bpse= sqrt(covmat[itar,itar])
		tp=c(tp,(bp/bpse)) 
		
	}
	
	if(tval>=0) pval = length(which(as.double(tp) >= as.double(tval)))/npermut
	if(tval<0) pval = length(which(as.double(tp) <  as.double(tval)))/npermut
	
	return(pval)
}

