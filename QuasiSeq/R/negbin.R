###########################################
### negbin was authored by Long Qu <long.qu@wright.edu> 
###########################################

### This extends the mgcv::negbin link function by providing different implementation of linkinv, initialize and an additional d2linkfun function 
negbin=function( link = "log", overdisp = stop("'overdisp' must be specified"))
{
	stopifnot(link=='log')
	
	ans=mgcv::negbin(1/overdisp,link)
	linkinv=function(eta) if(is.numeric(eta[1L])) pmin.int(pmax.int(exp(eta), .Machine$double.eps),.Machine$double.xmax)  else exp(eta) ## truncated exp function 
	environment(linkinv)=environment(ans$linkinv)
	ans$linkinv=ans$mu.eta=linkinv 

	dev.resids=function (y, mu, wt) 
	{
		Theta <- get(".Theta")[1]
		2 * wt * (y * log(pmax.int(1, y)/mu) - (y + Theta) * log((y + 
			Theta)/(mu + Theta)))
	}
	environment(dev.resids)=environment(ans$dev.resids)
	ans$dev.resids = dev.resids

	ans$d2linkfun=function(mu)-1/mu/mu ## used for observed info
	environment(ans$d2linkfun)=environment(ans$linkfun)
	
	ans$setTheta=function(theta)
	{
		assign('.Theta', theta, envir=parent.env(sys.frame(sys.nframe())))
	}
	environment(ans$setTheta)=environment(ans$getTheta)
	
	
	ans$initializers=list(
		expression({  ## lm + bias reduction on zeros only
			if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
			n <- rep(1, NROW(y))
			mustart = local({
				if(!is.null(offset)) offset=rep(0,  NROW(y))
				scale.fact=exp(offset-mean(offset))
				y=y/scale.fact
				if(sum(!duplicated(x))==NCOL(x)){
					  ans = pmax.int(0,lm.fit(x, y)$fitted.values) 
				}else ans = expm1(pmax.int(0,lm.fit(x, log1p(y))$fitted.values))
				y0=ans==0
				if(any(y0)) ans[y0]=ans[y0]+ .5/sum(y0)
				ans * scale.fact 
			})
		}),
	
	
		expression({ ## lm -> kmeans -> bias reduction
			if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
			n <- rep(1, NROW(y))
			mustart = local({
				if(!is.null(offset)) offset=rep(0,  NROW(y))
				scale.fact=exp(offset-mean(offset))
				y=y/scale.fact
				.tmp = lm.fit(x, log1p(y))
				.fit1=.tmp$fitted
				.km=kmeans(.fit1, sum(!is.na(coef(.tmp))))
				.n=.km$size[.km$cluster]
				.fitted = expm1(pmax.int(0,.km$centers)[.km$cluster])
				((1+.5/family$getTheta()/.n)*.fitted + .5/.n ) * scale.fact
			})
		}),

		expression(  ## this is pretty robust
		{  
			if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
			n=rep(1,NROW(y)); 
			.tmp=y+(y==0)/(2*sum(y==0)+1)
			mustart=.9*.tmp+.1*mean(.tmp)
		}),
	
		expression(  ## this is pretty robust
		{  
			if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
			n=rep(1,NROW(y)); 
			.tmp=y+(y==0)/(2*sum(y==0)+1)
			mustart=.5*.tmp+.5*mean(.tmp)
		})
	)
	attr(ans$initializers, 'i')=1L
	ans$initialize=ans$initializers[[1L]]
	
	if(!identical(environment(ans$variance), environment(ans$d3var))) environment(ans$d3var)=environment(ans$variance)  ## this seems to be a nearly harmless bug in mgcv_1.8-6 that prevents garbage collection
	
	class(ans)=c('fbrNBfamily', class(ans))
	ans
}

update.fbrNBfamily=function(object, overdisp, ...)
{  ## caution: this function has side effect of changing object
	tmp=object$setTheta(1/overdisp)
	object$family=sprintf("Negative Binomial(%.3f)", tmp)
	object
}
reinitialize.fbrNBfamily=function(object, ...)
{
	if(attr(object$initializers, 'i')==length(object$initializers)){
		warning('No more built-in initializers available')
		return(object)
	}
	attr(object$initializers, 'i')=attr(object$initializers, 'i')+1L
	object$initialize=object$initializers[[attr(object$initializers, 'i')]]
	object
}