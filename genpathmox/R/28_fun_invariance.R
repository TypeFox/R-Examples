#' @title Calculating the invariance test
#' @details
#' Internal function. \code{invariance} is called by \code{treemodel.pls}.
#' @param x matrix or data.frame with the data.
#' @param inner A square (lower triangular) boolean matrix representing the inner 
#' model (i.e. the path relationships between latent variables).
#' @param outer list of vectors with column indices or column names from Data indicating 
#' the sets of manifest variables forming each block (i.e. which manifest variables correspond to each block).
#' @param scaling optional list of string vectors indicating the type of 
#' measurement scale for each manifest variable specified in \code{blocks}.
#' \code{scaling} must be specified when working with non-metric variables.
#' Possible values: \code{"num"} (numeric), \code{"raw"}, \code{"nom"} (nominal), 
#' and \code{"ord"} (ordinal).
#' @param mode character vector indicating the type of measurement for each
#' block. Possible values are: \code{"A", "B", "newA", "PLScore", "PLScow"}. 
#' The length of \code{mode} must be equal to the length of \code{outer}.
#' @param scheme string indicating the type of inner weighting
#' scheme. Possible values are \code{"centroid"}, \code{"factorial"}, or
#' \code{"path"}.
#' @param scaled whether manifest variables should be standardized. 
#' Only used when \code{scaling = NULL}. When (\code{TRUE}, data is 
#' scaled to standardized values (mean=0 and variance=1). 
#' @return data frames containing the results of the invariance test 
#' @keywords internal
#' @export
	
invariance	<-	function(x,nodes,inner,outer,mode,scheme,scaling,scaled)
{
	lat 		= 	NULL
	block.h0	=	NULL
	blocks		=	NULL	
	
	for (i in 2 : length(nodes))
	{
		x.data	= 	x[nodes[[i]],]

		pls.node=plspm(x.data,inner,outer,mode,scaling,scheme,scaled=scaled)	


		x.node			=	pls.node$data
		inner.node		=	pls.node$model$IDM
		outer.node		=	pls.node$model$blocks
		
		
		latent.node = pls.node$scores


		
		lvs 			= 	nrow(inner.node)
		
		lat.node 		=	NULL
		block.node	=	NULL
		
		for(i in 1:ncol(inner.node)){lat.node = rbind(lat.node,as.matrix(latent.node[,i]))}
		
		if(scaled==TRUE)
		{
			for (k in 1:lvs) {block.node[[length(block.node)+1]] = as.matrix(cbind(1,scale(x.node[, outer.node[[k]]])))}
		}
		if(scaled==FALSE){
		for (k in 1:lvs) 
		{
			block.node[[length(block.node)+1]] = as.matrix(cbind(1,x.node[, outer.node[[k]]]))}
		} 	
 	
	 	block.h0 = rbind(block.h0,blockdiag(block.node)) 	
	 	
	 	blocks[[length(blocks)+1]]=blockdiag(block.node)
		lat = rbind(lat,lat.node)
	} 	
	
	block.h1 = blockdiag(blocks)
	
	lm0 = lm(lat~block.h0-1)
	lm1 = lm(lat~block.h1-1)
	
	SQR0 = sum(lm0$residuals^2)
	SQR1 = sum(lm1$residuals^2)
	dif.sqr = SQR0-SQR1
	
	df0 = (nrow(block.h0) - ncol(block.h0)) 
	df1 = (nrow(block.h1) - ncol(block.h1)) 
	dif.df = df0-df1
	p.val = pchisq(dif.sqr,dif.df,lower.tail=FALSE)  
	
	last = NULL
	for(i in 1:(length(outer)-1)){last[i] = tail(outer[[i]], n=1)+(i+1)}
	interpect = c(1,last)
	avg.weights = as.matrix(round(lm0$coefficients[-interpect],3))
	colnames(avg.weights) = "avg.weights"
	
	test = data.frame(chisq.statistic=dif.sqr,p.value=p.val,dfH0=df0,dfH1=df1)
	
	res = list(chisq.statistic=dif.sqr,
			p.value=p.val,
			dfH0=df0,
			dfH1=df1,
			avg.weights=avg.weights,
			test=test)
}

