#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.gect(): graph estimation via correlation thresholding (ct)       #                                                          #                           
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: Jul 15th 2011                                                   #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

##Main function
huge.ct = function(x, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL, verbose = TRUE)
{	
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
		x = scale(x)
		S = cor(x)
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
   				cat(paste(c("Conducting the graph estimation via correlation thresholding (ct) ....in progress:", floor(100*i/nlambda), "%"), collapse=""), "\r")
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
   				mes <- paste(c("Conducting the graph estimation via correlation thresholding (ct)....in progress:", floor(100*i/nlambda), "%"), collapse="")
   				cat(mes, "\r")
            	flush.console()
            }
		}
		fit$lambda = lambda
		rm(S,lambda)
		gc()
	}
		
	if(verbose)
	{
        cat("Conducting the graph estimation via correlation thresholding (ct)....done.             \r\n")
        flush.console()
    }
	gc()
	return(fit)
}