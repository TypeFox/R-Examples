print.mle <-
function(x,...)
{
	par.hat <- x$par.hat
	par.names <- names(unlist(formals(paste('r',x$dist,sep=''))))[-1]
	if(length(par.hat)!=length(par.names)) par.names <- par.names[1:length(par.hat)]
	par.se <- sqrt(diag(x$cov.hat))
	param <- rbind(Estimate=par.hat,Std.err=par.se)
	colnames(param) <- paste(par.names,'.hat',sep='')
	gof <- c(log.like=x$log.like,aic=x$aic,ad=x$ad,rho=x$rho)
	cat('-----------------------------------------\n')
	cat('       Maximum Likelihood Estimates\n')
	cat('-----------------------------------------\n')
	cat('Data object: ',x$data.name,'\n')
	cat('Distribution: ',x$dist,'\n')
	cat('\n')
	cat('--------- Parameter estimates -----------\n')
	cat('\n')
	print(param,digits=4)
	cat('\n')
	cat('---------- Goodness-of-Fit --------------\n')
	cat('\n')
	print(gof,digits=4)
	cat('-----------------------------------------\n')
}

