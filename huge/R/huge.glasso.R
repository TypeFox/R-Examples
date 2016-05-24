#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# glasso(): The graphical lasso (glasso) using sparse matrix output     #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: July 15th 2011                                                  #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

## Main function
huge.glasso = function(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, scr = NULL, cov.output = FALSE, verbose = TRUE){
	
	gcinfo(FALSE)
	n = nrow(x)
	d = ncol(x)
	fit = list()
	fit$cov.input = isSymmetric(x)
	if(fit$cov.input)
	{
		if(verbose) cat("The input is identified as the covriance matrix.\n")
		S = x
	}
	if(!fit$cov.input)
	{
		x = scale(x)
		S = cor(x)
	}
	rm(x)
	gc()
	if(is.null(scr)) scr = FALSE
	fit$scr = scr
	if(!is.null(lambda)) nlambda = length(lambda)
	if(is.null(lambda))
	{
		if(is.null(nlambda))
			nlambda = 10
		if(is.null(lambda.min.ratio))
			lambda.min.ratio = 0.1
		lambda.max = max(max(S-diag(d)),-min(S-diag(d)))
		lambda.min = lambda.min.ratio*lambda.max
		lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
	}
	
	fit$lambda = lambda
	fit$loglik = rep(-d,nlambda)
	fit$sparsity = rep(0,nlambda)
	fit$df = rep(0,nlambda)
	fit$path = list()
	fit$icov = list()
	fit$cov.output = cov.output	
	if(cov.output) fit$cov = list()
	out.glasso=NULL
	for(i in nlambda:1)
	{	
		z = which(rowSums(abs(S)>lambda[i])>1)
		q = length(z)

		if(q>0)
		{
			if(verbose){
				if(scr){
					cat(paste(c("Conducting the graphical lasso (glasso) wtih lossy screening....in progress:", floor(100*(1-i/nlambda)), "%"), collapse=""), "\r")
				}
				if(!scr){
					cat(paste(c("Conducting the graphical lasso (glasso) with lossless screening....in progress:", floor(100*(1-i/nlambda)), "%"), collapse=""), "\r")
				}		
				flush.console()
			}
			
			if(scr){
				if(!is.null(out.glasso))
					out.glasso = .C("hugeglassoscr", S = as.double(S[z,z]), W = as.double(tmp.cov[z,z]), T = as.double(tmp.icov[z,z]), dd = as.integer(q),lambda = as.double(lambda[i]), df = as.integer(0), PACKAGE="huge")
				if(is.null(out.glasso))
					out.glasso = .C("hugeglassoscr", S = as.double(S[z,z]), W = as.double(S[z,z]), T = as.double(diag(q)), dd = as.integer(q),lambda = as.double(lambda[i]), df = as.integer(0), PACKAGE="huge")
			}	
			
			else {
				if(!is.null(out.glasso))
					out.glasso = .C("hugeglasso", S = as.double(S[z,z]), W = as.double(tmp.cov[z,z]), T = as.double(tmp.icov[z,z]), dd = as.integer(q),lambda = as.double(lambda[i]), df = as.integer(0), PACKAGE="huge")
				if(is.null(out.glasso))
					out.glasso = .C("hugeglasso", S = as.double(S[z,z]), W = as.double(S[z,z]), T = as.double(diag(q)), dd = as.integer(q),lambda = as.double(lambda[i]), df = as.integer(0), PACKAGE="huge")
			}
			
			out.glasso$T = matrix(out.glasso$T,ncol=q)
			out.glasso$W = matrix(out.glasso$W,ncol=q)
		}
		if(q == 0) out.glasso = NULL

		tmp.icov = matrix(0,d,d);
		diag(tmp.icov) = 1/(diag(S)+lambda[i])
		tmp.cov = matrix(0,d,d);
		diag(tmp.cov) = diag(S)+lambda[i] 
		fit$path[[i]] = Matrix(0,d,d)
		if(!is.null(out.glasso))
		{
			tmp.icov[z,z] = out.glasso$T
			tmp.cov[z,z] = out.glasso$W
			fit$path[[i]][z,z] = abs(sign(out.glasso$T));
			diag(fit$path[[i]]) = 0;
			fit$sparsity[i] = as.double(out.glasso$df)/d/(d-1)
      fit$df[i] = out.glasso$df/2
			fit$loglik[i] = (log(det(out.glasso$T)) - sum(diag(out.glasso$T%*%S[z,z])) - (d-q))
		}
		fit$icov[[i]] = Matrix(tmp.icov)
		if(cov.output)
			fit$cov[[i]] = Matrix(tmp.cov)		
	}
	rm(S,out.glasso,tmp.cov,tmp.icov)
	gc()
	if(verbose){
   		cat("Conducting the graphical lasso (glasso)....done.                                          \r")
   		cat("\n")
        flush.console()
    }
    return(fit)
}