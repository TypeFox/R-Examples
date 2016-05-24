fhatnew <-
function(x,theta,mu,sigma){
	N=length(mu)
	n=length(x)
	dmult=matrix(0,N,n)
	for(j in 1:N){
	dmult[j,]=dnorm(x,mu[j],sigma)
	}
	f=theta%*%dmult
	mum=matrix(rep(mu,n),N,n,byrow=FALSE)
	xm=matrix(rep(x,N),N,n,byrow=TRUE)
	coeffg=(mum-xm)/sigma^2
	coeffgg=((mum-xm)^2/sigma^2-1)/sigma^2
	fg=theta%*%(dmult*coeffg)
	fgg=theta%*%(dmult*coeffgg)	
	return(cbind(t(f),t(fg),t(fgg)))
}
