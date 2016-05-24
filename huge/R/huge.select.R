#-------------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                   #
# huge.select(): Model selection using:                                   #
#			 	 1.rotation information criterion (ric)                   #
#                2.stability approach to regularization selection (stars) #
#                3.extended Bayesian informaition criterion (ebic)        #
# Authors: Tuo Zhao and Han Liu                                           #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                        #
# Date: Jan 21th 2012                                                    #
# Version: 1.2.0                                                          #
#-------------------------------------------------------------------------#

## Main Function
huge.select = function(est, criterion = NULL, ebic.gamma = 0.5, stars.thresh = 0.1, stars.subsample.ratio = NULL, rep.num = 20, verbose = TRUE){

	gcinfo(FALSE)
	
	if(est$cov.input){
		cat("Model selection is not available when using the covariance matrix as input.")
		class(est) = "select"
    	return(est)
	}
	if(!est$cov.input)
	{		
		if(est$method == "mb"&&is.null(criterion))
			criterion = "ric"
		if(est$method == "ct"&&is.null(criterion))
			criterion = "stars"
		if(est$method == "glasso"&&is.null(criterion))
			criterion = "ebic"
	
		n = nrow(est$data)
		d = ncol(est$data)
		nlambda = length(est$lambda)	
	
		if(criterion == "ric")
		{
			if(verbose)
			{
				cat("Conducting rotation information criterion (ric) selection....")
        		flush.console()
        	}
		
			if(n>rep.num){
				nr = rep.num
				r = sample(n,rep.num)
			}
			if(n<=rep.num){
				nr = n
				r = 1:n
			} 
			
			out=.C("RIC",X = as.double(est$data),dd = as.integer(d),nn=as.integer(n),r=as.integer(r),nr=as.integer(nr),lambda_opt = as.double(0),PACKAGE="huge")
			est$opt.lambda = out$lambda_opt/n
			rm(out)
			gc()
		
			if(verbose){
				cat("done\n")
				flush.console()
			}
		
			if(verbose)
			{
				cat("Computing the optimal graph....")
				flush.console()
			}
			
			if(est$opt.lambda>max(cor(est$data)))
				est$refit = Matrix(0,d,d)
			else{
		
			if(est$method == "mb")
				est$refit = huge.mb(est$data, lambda = est$opt.lambda, sym = est$sym, idx.mat = est$idx.mat, verbose = FALSE)$path[[1]]
			if(est$method == "glasso")
			{
				if(!is.null(est$cov))
				{
					tmp = huge.glasso(est$data, lambda = est$opt.lambda, scr = est$scr, cov.output = TRUE, verbose = FALSE)
					est$opt.cov = tmp$cov[[1]]
				}
				if(is.null(est$cov))
					tmp = huge.glasso(est$data, lambda = est$opt.lambda, verbose = FALSE)
			
				est$refit = tmp$path[[1]]
				est$opt.icov = tmp$icov[[1]]
				rm(tmp)
				gc()
			}
			if(est$method == "ct")
				est$refit = huge.ct(est$data, lambda = est$opt.lambda, verbose = FALSE)$path[[1]]
			}
			est$opt.sparsity=sum(est$refit)/d/(d-1)
		
			if(verbose){
				cat("done\n")
				flush.console()
			}
		}
	
		if(criterion == "ebic"&&est$method == "glasso")
		{
			if(verbose)
			{
				cat("Conducting extended Bayesian information criterion (ebic) selection....")
        		flush.console()
      }
			est$ebic.score = -n*est$loglik + log(n)*est$df + 4*ebic.gamma*log(d)*est$df
			est$opt.index = which.min(est$ebic.score)
			est$refit = est$path[[est$opt.index]]
			est$opt.icov = est$icov[[est$opt.index]]
			if(est$cov.output)
				est$opt.cov = est$cov[[est$opt.index]]
  			est$opt.lambda = est$lambda[est$opt.index]
  			est$opt.sparsity = est$sparsity[est$opt.index]
  			if(verbose){
				cat("done\n")
				flush.console()
			}
		}		
	
		if(criterion == "stars"){
			if(is.null(stars.subsample.ratio))
			{
				if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
				if(n<=144) stars.subsample.ratio = 0.8
			} 
	
			est$merge = list()
			for(i in 1:nlambda) est$merge[[i]] = Matrix(0,d,d)
		
  			for(i in 1:rep.num)
  			{
  				if(verbose)
  				{
					mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/rep.num), "%"), collapse="")
   					cat(mes, "\r")
            		flush.console()	
				}
    			ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
    		
    			if(est$method == "mb")
    				tmp = huge.mb(est$data[ind.sample,],lambda = est$lambda, scr = est$scr, idx.mat = est$idx.mat, sym = est$sym, verbose = FALSE)$path
       			if(est$method == "ct")
    				tmp = huge.ct(est$data[ind.sample,], lambda = est$lambda,verbose = FALSE)$path
    			if(est$method == "glasso")
    				tmp = huge.glasso(est$data[ind.sample,], lambda = est$lambda, scr = est$scr, verbose = FALSE)$path
    			
    			for(i in 1:nlambda)
    				est$merge[[i]] = est$merge[[i]] + tmp[[i]]
    		
    			rm(ind.sample,tmp)
   				gc()
			}
		
			if(verbose){
				mes = "Conducting Subsampling....done.                 "
      		cat(mes, "\r")
      		cat("\n")
      		flush.console()
    		}
        
    		est$variability = rep(0,nlambda)
			for(i in 1:nlambda){
				est$merge[[i]] = est$merge[[i]]/rep.num
    			est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(d*(d-1))
    		}
    
    		est$opt.index = max(which.max(est$variability >= stars.thresh)[1]-1,1)
   			est$refit = est$path[[est$opt.index]]
  			est$opt.lambda = est$lambda[est$opt.index]
  			est$opt.sparsity = est$sparsity[est$opt.index]
  			if(est$method == "glasso")
  			{
  				est$opt.icov = est$icov[[est$opt.index]]
				if(!is.null(est$cov))
					est$opt.cov = est$cov[[est$opt.index]]
  			}
		}
    	est$criterion = criterion
  		class(est) = "select"
    	return(est)
    }  	
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

print.select = function(x, ...)
{
	if(x$cov.input){
		cat("Model selection is not available when using the covariance matrix as input.")
	}
	if(!x$cov.input)
	{
		if(x$method == "ct")
			cat("Model: graph gstimation via correlation thresholding(ct)\n")
		if(x$method == "glasso")
			cat("Model: graphical lasso (glasso)\n")
		if(x$method == "mb")
			cat("Model: Meinshausen & Buhlmann Graph Estimation (mb)\n")

		cat("selection criterion:",x$criterion,"\n")
		if((x$method != "ct")&&x$scr)
			cat("lossy screening: on\n")
		cat("Graph dimension:",ncol(x$data),"\n")
		cat("sparsity level", x$opt.sparsity,"\n")
	}
}

plot.select = function(x, ...){
	if(x$cov.input){
		cat("Model selection is not available when using the covariance matrix as input.")
	}
	if(!x$cov.input)
	{
		par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    
    	g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
		layout.grid = layout.fruchterman.reingold(g)
	
		plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA)	  
    	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
    	lines(x$opt.lambda,x$opt.sparsity,type = "p")  
    }
}