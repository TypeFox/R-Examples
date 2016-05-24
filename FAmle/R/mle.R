mle <-
function(x,dist,start=NULL,method='Nelder-Mead')
{
	if(is.null(start))
		stop('\'start\' must be specified!')
	suppressWarnings(fit <- optim(start,function(k) -sum(distr(x,dist,k,'d',log=TRUE)),hessian=TRUE,method=method))
	par.hat <- fit$par
	cov.hat <- solve(fit$hessian)
	if(sum(is.na(suppressWarnings(sqrt(diag(cov.hat)))))!=0)
		stop('Try with other starting values!')
	n <- length(x)
	i <- 1:n
	x.info <- cbind(i=i,x=x,z=sort(x),Fx=distr(x,dist,par.hat,'p'),Fz=distr(sort(x),dist,par.hat,'p'),
		Emp=i/(n+1),zF=distr(i/(n+1),dist,par.hat,'q'),fx=distr(x,dist,par.hat,'d'),
		fz=distr(sort(x),dist,par.hat,'d'))
	aic <- 2*length(par.hat) + 2*fit$value
	rho <- cor(sort(x),x.info[,'zF'])
	Fz <- x.info[,'Fz']
	ad <- -n-sum((2*i-1)/n*(log(Fz)+log(1-Fz[n+1-i])))
	out <- list(fit=fit,x.info=x.info,dist=dist,par.hat=par.hat,cov.hat=cov.hat,k=length(par.hat),
		n=n,log.like=-fit$value,aic=aic,ad=ad,data.name=deparse(substitute(x)),rho=rho)
	class(out) <- 'mle'
	return(out)
}
