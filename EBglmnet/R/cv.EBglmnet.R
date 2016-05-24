cv.EBglmnet <-function(x,y,family=c("gaussian","binomial"),prior= c("lassoNEG","lasso","elastic net"),nfolds=5,foldId,Epis = FALSE,group = FALSE, verbose = 0 ){
verbose = verbose -1;
	if(prior!="lassoNEG" && group)
	{
		warning("group EBlasso is only available on the NEG prior, ignore this parameter")
	}
	if(!Epis && group)
	{
		warning("group EBlasso is designed for epistasis analysis, ignore this parameter")
		group = FALSE;
	}
	if (missing(foldId)) 
	{
		N = nrow(x);
		if(N%%nfolds!=0){
			foldId 			= sample(c(rep(1:nfolds,floor(N/nfolds)),1:(N%%nfolds)),N);
		}else{
			foldId 			= sample(rep(1:nfolds,floor(N/nfolds)),N);
		}
	}
	this.call=match.call()#returns a call in which all of the specified arguments are specified by their full names.
	family=match.arg(family);
	prior = match.arg(prior);
	#check parameters
	y=drop(y) # we dont like matrix responses unless we need them
	np=dim(x)
	if(is.null(np)|(np[2]<=1))stop("x should be a matrix with 2 or more columns")
	nobs=as.integer(np[1])
	dimy=dim(y)
	nrowy=ifelse(is.null(dimy),length(y),dimy[1])
	if(nrowy!=nobs)stop(paste("number of observations in y (",nrowy,") not equal to the number of rows of x (",nobs,")",sep=""))
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
		cv=switch(family,
		"gaussian"=EBelasticNet.GaussianCV(x,y,nfolds,foldId,Epis),
		"binomial"=EBelasticNet.BinomialCV(x,y,nfolds,foldId,Epis)
		)
		opt_para = cv$optimal;
		alpha = opt_para[1];
		lambda = opt_para[2];
		
		fit=switch(family,
		"gaussian"=EBelasticNet.Gaussian(x,y,lambda,alpha,Epis,verbose),
		"binomial"=EBelasticNet.Binomial(x,y,lambda,alpha,Epis,verbose)
		)

		
	}else if(prior=="lasso")
	{
		cv=switch(family,
		"gaussian"=EBlassoNE.GaussianCV(x,y,nfolds,foldId,Epis,verbose),
		"binomial"=EBlassoNE.BinomialCV(x,y,nfolds,foldId,Epis,verbose)
		)
		alpha 	= 1;
		lambda 	= cv$optimal
		fit=switch(family,
		"gaussian"=EBelasticNet.Gaussian(x,y,lambda,alpha,Epis,verbose),
		"binomial"=EBelasticNet.Binomial(x,y,lambda,alpha,Epis,verbose)
		)
		
	}else
	{
		cv=switch(family,
		"gaussian"=EBlassoNEG.GaussianCV(x,y,nfolds,foldId, Epis,verbose,group),
		"binomial"=EBlassoNEG.BinomialCV(x,y,nfolds,foldId, Epis,verbose,group)
		)	
		#EBlassoNEG.BinomialCV <-function(BASIS,Target,nFolds,foldId,Epis=FALSE,verbose = 0, group = FALSE)
		opt_para = cv$optimal;
		a = opt_para[1];
		b = opt_para[2];
		fit=switch(family,
		"gaussian"=EBlassoNEG.Gaussian(x,y,a,b,Epis,verbose,group),
		"binomial"=EBlassoNEG.Binomial(x,y,a,b,Epis,verbose,group)
		)
	}
	#output = list();
	#output$CrossValidation = cv$CrossValidation;
	#output$fit = fit;
	output = c(cv,fit);
	output$family = family
	output$prior = prior	
	output$call=this.call
	output$nobs=nobs
	output$nfolds = nfolds
	class(output)="cv.EBglmnet"
	return(output)	
}
