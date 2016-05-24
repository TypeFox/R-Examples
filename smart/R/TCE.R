TCE <-
function(x, method, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL, verbose = TRUE)
{	
  if(!method%in%c("pearson","kendall","spearman","npn","ns")){ 
      cat("method not provided correctly, should be one of (pearson,kendall,spearman,npn,ns) \n")
      return(NULL)
  }

	gcinfo(FALSE)
	n = nrow(x);
	d = ncol(x);
	fit = list()
	fit$cov.input = isSymmetric(x);
	if(fit$cov.input)
	{
		if(verbose) cat("The input is identified as the covriance matrix.\n")
		S = cov2cor(x);
	}
	if(!fit$cov.input)
	{
     tmp = smart.npn(x,npn.func=method,verbose=verbose)
     S = tmp$cov
	}
	
	rm(x)
	gc()
	diag(S) = 0
	S = abs(S)
  	S.rank = order(S,decreasing = TRUE)
	gc()
 		
	if(is.null(lambda))
	{
		if(is.null(nlambda))
			nlambda = 20
		if(is.null(lambda.min.ratio))
			lambda.min.ratio = 0.05
		 		
 		density.max = lambda.min.ratio*d*(d-1)/2
 		density.min = 1
 		density.all = ceiling(seq(density.min,density.max,length = nlambda))*2
 		fit$sparsity = density.all/d/(d-1)
 		fit$lambda = S[S.rank[density.all]]
 		rm(density.max,lambda.min.ratio,density.min,S)
		gc()
 						
 		fit$path = list()
		for(i in 1:nlambda)
		{
			fit$path[[i]] = Matrix(0,d,d)
			fit$path[[i]][S.rank[1:density.all[i]]] = 1
			if(verbose)
			{
   				cat(paste(c("Conducting the transelliptical correlation estimation:", floor(100*i/nlambda), "%"), collapse=""), "\r")
            	flush.console()
            }	
		}
		rm(density.all,nlambda,S.rank)
		gc()
	}

	if(!is.null(lambda))
	{
		nlambda = length(lambda)
		fit$path = list()
		fit$sparsity = rep(0,nlambda)
		for(i in 1:nlambda)
		{
			fit$path[[i]] = Matrix(0,d,d)
			fit$path[[i]][S > lambda[i]] = 1
			fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
			if(verbose)
			{
   				mes <- paste(c("Conducting the transelliptical correlation estimation:", floor(100*i/nlambda), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()
            }
		}
		fit$lambda = lambda
		rm(S,lambda)
		gc()
	}
	fit$method = method	
		
	if(verbose)
	{
        cat("done \n")
        flush.console()
    }
	gc()
	
  class(fit) = "TCE"
	return(fit)
}
