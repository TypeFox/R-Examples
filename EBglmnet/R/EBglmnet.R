EBglmnet <-function(x,y,family=c("gaussian","binomial"),prior= c("lassoNEG","lasso","elastic net"),hyperparameters,Epis = FALSE,group = FALSE, verbose = 0 ){
	if(prior!="lassoNEG" && group)
	{
		warning("group EBlasso is only available on the NEG prior, ignore this parameter")
	}
	if(!Epis && group)
	{
		warning("group EBlasso is designed for epistasis analysis, ignore this parameter")
		group = FALSE;
	}
	if(missing(hyperparameters))
	{
		warning("hyperparameters controling number of nonzero effects need to be specified; run cv.EBglmnet to determine the parameters \n")
		hyperparameters = c(0.5, 0.5);
	}
	
	this.call=match.call()#returns a call in which all of the specified arguments are specified by their full names.
	family=match.arg(family);
	prior = match.arg(prior); #1) no prior specified: default the 1st one; 2) mismatched prior names: error reported.
	y=drop(y) # we dont like matrix responses unless we need them
	np=dim(x)
	if(is.null(np)|(np[2]<=1))stop("x should be a matrix with 2 or more columns")
	nobs=as.integer(np[1])
	dimy=dim(y)
	nrowy=ifelse(is.null(dimy),length(y),dimy[1])
	if(nrowy!=nobs)stop(paste("number of observations in y (",nrowy,") not equal to the number of rows of x (",nobs,")",sep=""))
  	#
	if(family=="binomial")
	{
	    ## Need to construct a y matrix, and include the weights
		y=as.factor(y)
		ntab=table(y)
		minclass=min(ntab)
		if(minclass<=1)stop("one binomial class has 1 or 0 observations; not allowed")
		if(minclass<8)warning("one binomial class has fewer than 8  observations; dangerous ground")

		nc=as.integer(length(ntab))
		y0=diag(nc)[as.numeric(y),]
		y = y0[,2];
	}
	#end check
	if(prior=="elastic net")
	{
		alpha = hyperparameters[1];
		lambda = hyperparameters[2];
		if(alpha>1 || alpha<0)
		{
			warning("EBEN: alpha not in range of [0,1]; set to 1")
			alpha=1
		}
		if(lambda<0)
		{
			warning("EBEN: lambda should be a positive number; set lambda = 1")
			lambda = 0.1;
		}		
		
		fit=switch(family,
		"gaussian"=EBelasticNet.Gaussian(x,y,lambda,alpha,Epis,verbose),
		"binomial"=EBelasticNet.Binomial(x,y,lambda,alpha,Epis,verbose)
		)		
	}else if(prior=="lasso")
	{
		alpha = 1;
		lambda = hyperparameters[1];
		if(lambda<0)
		{
			warning("EBEN: lambda should be a positive number; set lambda = 1")
			lambda = 0.1;
		}
		
		fit=switch(family,
		"gaussian"=EBelasticNet.Gaussian(x,y,lambda,alpha,Epis,verbose),
		"binomial"=EBelasticNet.Binomial(x,y,lambda,alpha,Epis,verbose)
		)		
	}else
	{
		a = hyperparameters[1];
		b = hyperparameters[2];
		if(a<=-1.5)
		{
			warning("EBlassoNEG has support of a> 1.5 and b>0, set a = 0.1")
			a = 0.1;
		}	
		if(b<=0)
		{
			warning("EBlassoNEG has support of a> 1.5 and b>0, set b = 0.1")
			b = 0.1;
		}
		
		fit=switch(family,
		"gaussian"=EBlassoNEG.Gaussian(x,y,a,b,Epis,verbose,group),
		"binomial"=EBlassoNEG.Binomial(x,y,a,b,Epis,verbose,group)
		)		
	}
	fit$family = family
	fit$prior = prior
	fit$call=this.call
	fit$nobs=nobs
	class(fit)=c(class(fit),"glmnet")
	return(fit)	
}
