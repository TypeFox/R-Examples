#-----------------------------------------------------------------------#
# Package: High-dimensional Undirecte d Graph Estimation                #
# huge(): The user interface for huge.mb(), huge.glasso() and huge.ct() #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: Jul 15th 2011                                                   #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

## Main function
huge = function(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, method = "mb", scr = NULL, scr.num = NULL, cov.output = FALSE, sym = "or", verbose = TRUE)
{	
	gcinfo(FALSE)
	est = list()
	est$method = method	
		
	if(method == "ct")
	{
		fit = huge.ct(x, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$cov.input = fit$cov.input
		rm(fit)
		gc()
	}
	
	if(method == "mb")
	{	
		fit = huge.mb(x, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, scr = scr, scr.num = scr.num, sym = sym, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$df = fit$df
		est$idx_mat = fit$idx_mat
		est$sym = sym
		est$scr = fit$scr
		est$cov.input = fit$cov.input
		rm(fit,sym)
		gc()
	}
	
	
	if(method == "glasso")
	{
		fit = huge.glasso(x, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, scr = scr, cov.output = cov.output, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$icov = fit$icov
		est$df = fit$df
		est$sparsity = fit$sparsity
		est$loglik = fit$loglik
		if(cov.output)
			est$cov = fit$cov
		est$cov.input = fit$cov.input
		est$cov.output = fit$cov.output	
		est$scr = fit$scr
		rm(fit)
		gc()
	}			
			
	est$data = x
	
	rm(x,scr,lambda,lambda.min.ratio,nlambda,cov.output,verbose)
	gc()
	class(est) = "huge"
	return(est)
}

print.huge = function(x, ...)
{	
	if(x$method == "ct")
		cat("Model: graph estimation via correlation thresholding (ct)\n")
	if(x$method == "glasso")
		cat("Model: graphical lasso (glasso)\n")
	if(x$method == "mb")
		cat("Model: Meinshausen & Buhlmann graph estimation (mb)\n")
	
	if((x$method != "ct")&&(x$scr)) cat("lossy screening: on\n")

	if(x$cov.input) cat("Input: The Covariance Matrix\n")
	if(!x$cov.input) cat("Input: The Data Matrix\n")
	
	cat("Path length:",length(x$lambda),"\n")
	cat("Graph dimension:",ncol(x$data),"\n")
	cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}

plot.huge = function(x, align = FALSE, ...){
  gcinfo(FALSE)
	
	if(length(x$lambda) == 1)	par(mfrow = c(1, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$lambda) == 2)	par(mfrow = c(1, 3), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$lambda) >= 3)	par(mfrow = c(1, 4), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	
	if(length(x$lambda) <= 3)	z.final = 1:length(x$lambda)
	
	if(length(x$lambda) >=4){
		z.max = max(x$sparsity)
		z.min = min(x$sparsity)
		z = z.max - z.min
		z.unique = unique(c(which(x$sparsity>=(z.min + 0.03*z))[1],which(x$sparsity>=(z.min + 0.07*z))[1],which(x$sparsity>=(z.min + 0.15*z))[1]))

		
		if(length(z.unique) == 1){
			if(z.unique<(length(x$lambda)-1))	z.final = c(z.unique,z.unique+1,z.unique+2)
			if(z.unique==(length(x$lambda)-1)) z.final = c(z.unique-1,z.unique,z.unique+1)
			if(z.unique==length(x$lambda)) 	z.final = c(z.unique-2,z.unique-1,z.unique)
		}
		
		if(length(z.unique) == 2){
			if(diff(z.unique)==1){
				if(z.unique[2]<length(x$lambda)) z.final = c(z.unique,z.unique[2]+1) 
				if(z.unique[2]==length(x$lambda)) z.final = c(z.unique[1]-1,z.unique)
			}
			if(diff(z.unique)>1) z.final = c(z.unique[1],z.unique[1]+1,z.unique[2])
		}
		
		if(length(z.unique) == 3) z.final = z.unique
		
		rm(z.max,z.min,z,z.unique)
		gc()
		
	}
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
	
	lines(x$lambda[z.final],x$sparsity[z.final],type = "p")
	
	if(align){
		layout.grid = layout.fruchterman.reingold(graph.adjacency(as.matrix(x$path[[z.final[length(z.final)]]]), mode="undirected", diag=FALSE))
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g)
			gc()
		}
	rm(layout.grid)
	}
	if(!align){
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			layout.grid = layout.fruchterman.reingold(g)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g,layout.grid)
			gc()
		}
	}
	if(align) cat("Three plotted graphs are aligned according to the third graph\n")
}