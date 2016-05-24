#'@title Invariance Test 
#'
#'@description
#'The invariance test is a test that allows to verify the existence of
#'common weights for the different local PLS-PM models identified by one or more 
#'segmentation variable.
#'
#'@details
#' The \code{"x"} refers to a matrix or a data.farme that contains all individuals 
#' used for the global PLS-PM estimation
#' The \code{"nodes"} is a list of vectors. Each vector contains the position of the  
#' individual that belongs to a specific node. The position is idenfied by the number of row.
#' For example, the row 4 corresponds to the individual 4 
#' The other parameters are the classical parameter of the fucntion \code{"plspm"}
#'
#' @param x Matrix or data frame containing the manifest variables.
#' @param nodes List of vectors. Each vector contains the position of the  
#' individual that belongs to a specific node.
#' @param inner A square (lower triangular) boolean matrix representing 
#' the inner model (i.e. the path relationships between latent variables).
#' @param outer list of vectors with column indices or column names
#' from \code{Data} indicating the sets of manifest variables forming 
#' each block (i.e. which manifest variables correspond to each block).
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
#' 
#' @return An data.farme \code{res}. Basically a list with the
#' following results:
#' @return \item{chisq.statistic}{A Number; X^2 statistic}
#' @return \item{p.value}{A Number; p-value}
#' @return \item{dfH0}{A Number; degree of freedom null Hypotheis}
#' @return \item{dfH1}{A Number; degree of freedom alternative Hypotheis}
#' @return \item{avg.weights}{data frame of the common weights if they exist}
#' @return \item{test}{data frame with summry information of the invariance test}
#' 
#'
#' @author Giuseppe Lamberti, Tomas Aluja
#' 
#' 
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation. 
#'
#' @examples
#'  \dontrun{
#'  ## example of PLS-PM in alumni satisfaction
#'  
#'  data(fibtele)
#'  
#'  # select manifest variables
#'  data.fib <-fibtele[,12:35]
#'  
#'  # define inner model matrix
#'  Image 			= rep(0,5)
#'	 Qual.spec	= rep(0,5)
#'	 Qual.gen		= rep(0,5)
#'	 Value			= c(1,1,1,0,0)
#'	 Satis			= c(1,1,1,1,0)
#'  inner.fib <- rbind(Image,Qual.spec, Qual.gen, Value, Satis)
#'  colnames(inner.fib) <- rownames(inner.fib)
#'  
#'  # blocks of indicators (outer model)
#'  outer.fib <- list(1:8,9:11,12:16,17:20,21:24)
#'  modes.fib  = rep("A", 5)
#'  
#'  # apply plspm
#'  pls.fib <- plspm(data.fib, inner.fib, outer.fib, modes.fib)
#'                  
#'
#'  # re-ordering those segmentation variables with ordinal scale 
#'   seg.fib= fibtele[,2:11]
#'  
#'	seg.fib$Age = factor(seg.fib$Age, ordered=T)
#'	seg.fib$Salary = factor(seg.fib$Salary, 
#'		levels=c("<18k","25k","35k","45k",">45k"), ordered=T)
#'	seg.fib$Accgrade = factor(seg.fib$Accgrade, 
#'		levels=c("accnote<7","7-8accnote","accnote>8"), ordered=T)
#'	seg.fib$Grade = factor(seg.fib$Grade, 
#'		levels=c("<6.5note","6.5-7note","7-7.5note",">7.5note"), ordered=T)
#'
#'  # Pathmox Analysis
#'  fib.pathmox=pls.pathmox(pls.fib,seg.fib,signif=0.05,deep=2,size=0.2,n.node=20)
#'
#'  # Select the terminal nodes
#'  ls(fib.pathmox)
#' 
#'  terminal.nodes=fib.pathmox$terminal 
#'
#'  # Invariance test
#'  inv.test=invariance_test(data.fib,terminal.nodes,inner.fib,
#'				outer.fib,modes.fib,scheme="centroid",scale=T)
#'  inv.test
#'  }
#'
invariance_test	<-	function(x,nodes,inner,outer,mode,scheme,scaling,scaled)
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
		
		if (scale == FALSE){latent.node = pls.node$scores}
		if (scale == TRUE){latent.node = pls.node$latents}

		
		lvs 			= 	nrow(inner.node)
		
		lat.node 		=	NULL
		block.node	=	NULL
		
		for(i in 1:ncol(inner.node)){lat.node = rbind(lat.node,as.matrix(latent.node[,i]))}
		for (k in 1:lvs) {block.node[[length(block.node)+1]] = as.matrix(cbind(1,x.node[, outer.node[[k]]]))}
	 	
 	
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

