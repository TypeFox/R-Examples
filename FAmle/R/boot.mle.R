boot.mle <-
function(model,B=200,seed=NULL,start=NULL,method='Nelder-Mead')
{
	if(!is.null(seed)) set.seed(seed)
	if(is.null(start)) start <- model$par.hat
	param <- array(NaN,c(B,model$k))
	gof <- array(NaN,c(B,2),dimnames=list(NULL,c('ad','rho')))
	t1 <- Sys.time()
	for(i in 1:B)
	{
		fit <- NULL
		x <- distr(x=model$n,model=model,type='r')
		try(fit <- mle(x,model$dist,start,method=method),silent=TRUE)
		if(!is.null(fit))
		{
			param[i,] <- fit$par.hat
			gof[i,] <- c(fit$ad,fit$rho)
		}
	}
	tot.time <- Sys.time() - t1
	colnames(param) <- paste(names(formals(paste('r',model$dist,sep='')))[2:(1+model$k)],'.star',sep='')
	p.value <- c(ad=sum(gof[,'ad']>=model$ad,na.rm=TRUE),rho=sum(gof[,'rho']<=model$rho,na.rm=TRUE))/(B+1)
	failure.rate <- mean(is.na(param[,1]))
	out <- list(model=model,B=B,seed=seed,par.star=param,gof=gof,
		p.value=p.value,failure.rate=failure.rate,total.time=tot.time)
	return(out)
}