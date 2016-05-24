

metropolis <-
function(model,iter=1000,tun=2,trans.list=NULL,start=NULL,variance=NULL,prior=NULL,burn=0,uniroot.interval=c(-100,100),
	pass.down.to.C=FALSE)
{
	if(class(model)!='mle' & (is.null(start) || is.null(variance)))
		stop('If argument \'model\' is not from class \'mle\', both \'start\' and \'variance\' must be specified!')
	if(class(model)=='mle')
	{
		x <- model$x.info[,'x']
		k <- model$k
	}
	else
	{
		x <- model$x
		k <- length(start)
		if(!is.null(trans.list))
		{
			model$par.hat <- vector()
			for(i in 1:k) model$par.hat[i] <- trans.list[[i]](start[i])
		}
		else model$par.hat <- start
	}
	if(!is.null(trans.list) & length(trans.list) != k)
		stop('the length of \'trans.list\' must match the number of unknown parameters!')
	else if(is.null(trans.list)) trans.list <- lapply(as.list(1:k),function(g) function(x) x)
	log.like <- function(param)
	{
		param <- as.vector(param)
		for(i in 1:k) param[i] <- trans.list[[i]](param[i])
		distr(x,model$dist,param,'d',log=TRUE)
	}
	prior.yes <- 'yes'
	if(is.null(prior))
	{
		prior <- function(x) rep(1,length(x))
		prior.yes <- 'no'
	}
	fit <- model$fit
	if(!is.null(trans.list))
	{
		if(class(model)=='mle')
		{
			start.trans <- sapply(as.list(1:k),function(h) uniroot(function(g)
				trans.list[[h]](g)-model$par.hat[h],uniroot.interval)$root)
			fit <- optim(start.trans,function(g) -sum(log.like(g)),hessian=TRUE)
		}
		else 	start.trans <- sapply(as.list(1:k),function(h) uniroot(function(g)
			trans.list[[h]](g)-model$par.hat[h],uniroot.interval)$root)
	
	}
	if(is.null(start)) M <- fit$par
	else M <- start
	if(is.null(variance)) V <- solve(fit$hessian)*tun
	else V <- variance*tun
	ratio.test <- function(a,b) exp(sum(log.like(b)-log.like(a))+sum(log(prior(b))-log(prior(a))))
	t1 <- Sys.time()
	if(pass.down.to.C)
	{
		out.C <- .Call('metropolis',as.integer(iter),as.double(M),V,quote(ratio.test(a,b)),new.env(),PACKAGE='FAmle')
		sims <- t(out.C[[1]])
		rate <- out.C[[2]]
	}
	else
	{
		sims <- array(NaN,c(iter,k))
		sims[1,] <- M
		rate <- rep(0,iter)
		for(i in 2:iter)
		{
			sims[i,] <- rmvnorm(1,sims[i-1,],V)
			suppressWarnings(ratio <- ratio.test(sims[i-1,],sims[i,]))
			if(!is.na(ratio) & runif(1) < ratio) rate[i] <- 1
			else sims[i,] <- sims[i-1,]
		}
	}
	t2 <- difftime(Sys.time(),t1,units='mins')
	sims.out <- sapply(as.list(1:k),function(g) trans.list[[g]](sims[,g]))
	colnames(sims.out) <- names(formals(paste('r',model$dist,sep='')))[2:(1+k)]
	if(burn!=0) sims.burnt <- sims.out[-c(1:burn),]
	else sims.burnt <- sims.out
	out <- list(rate=mean(rate),total.time=t2,sims.all=sims.out,sims=sims.burnt,
		input=model,iter=iter,prior=prior.yes,burn=burn,M=sims[iter,],V=cov(sims[-c(1:burn),]))
	class(out) <- 'metropolis'
	return(out)
}
