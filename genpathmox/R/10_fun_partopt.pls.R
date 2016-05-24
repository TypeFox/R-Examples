
#' @title Defining the optimum partition given a set of segmentation variables
#' @details
#' Internal function. \code{partopt.pls} is called by \code{pls.pathmox}.
#' @param x matrix or data frame containing the data.
#' @param y matrix or data.frame of the segmentation variables.
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
#' @param method string indicating the method: LM or LAD
#' @param n.node number indicating a stop condition
#' @param \dots Further arguments passed on to \code{\link{partopt.pls}}. 
#' @return list containing information of the optimum partition given a set of segmentation variables
#' @keywords internal
#' @export

partopt.pls	<-	function(x,y,inner,outer,mode,scheme,scaling,scaled,method,n.node,...)
{
	a.p	=	all.part.pls(x,y,inner,outer,mode,scheme,scaling,scaled,method,n.node)
	
	if (any(!is.null(a.p$pvl))){
		
	fun.min	=	f.min(a.p$pvl)
	
		p.bin.opt		=	a.p$p.bin[[fun.min$p.min]]                  
		level.opt		=	a.p$level[[fun.min$p.min]] 
		variable.opt	=	a.p$variable[[fun.min$p.min]]  
		pvl.opt		=	a.p$pvl[fun.min$p.min]
		Ftest.opt		=	a.p$Ftest[fun.min$p.min]       
		g1.ind.opt	=	a.p$g1.ind[fun.min$p.min]      
		g2.ind.opt	=	a.p$g2.ind[fun.min$p.min]       
		modtwo.opt	=	a.p$modtwo[[fun.min$p.min]]		
		df0.opt		=	a.p$df0[fun.min$p.min]
		df1.opt		=	a.p$df1[fun.min$p.min]
		
		variable		=	a.p$variable[fun.min$all.v]   
		Ftest		=	a.p$Ftest[!is.na(a.p$Ftest)]
		pvl			=	a.p$pvl[!is.na(a.p$pvl)]
		df0			=	a.p$df0[!is.na(a.p$df0)]
		df1			=	a.p$df1[!is.na(a.p$df1)]
		g1.ind		=	a.p$g1.ind[!is.na(a.p$g1.ind)]
		g2.ind		=	a.p$g2.ind[!is.na(a.p$g2.ind)]
		
		i = 1
		level = NULL
		
		while (i < length(unlist(a.p$level)))
		{
			level = rbind(level,unlist(a.p$level)[c(i,i+1)])
			i = i+2
		}
		colnames(level)	=	c("g1.mod","g2.mod")

		candidates	=	data.frame(variable,Ftest,pvl,df0,df1,g1.ind,g2.ind,level)
		candidates	=	candidates[order(candidates[,2],decreasing=T),]

	list(candidates=candidates,
		p.bin.opt=p.bin.opt,
		variable.opt=variable.opt,
		level.opt=level.opt,
		Ftest.opt=Ftest.opt,
		pvl.opt=pvl.opt,
		indg1.opt=g1.ind.opt,
		indg2.opt=g2.ind.opt,
		df0.opt=df0.opt,
		df1.opt=df1.opt,
		modtwo.opt=modtwo.opt)
	}
	else
	{
		list(candidates=NULL,
		p.bin.opt=NULL,
		variable.opt=NULL,
		level.opt=NULL,
		Ftest.opt=NULL,
		pvl.opt=NULL,
		indg1.opt=NULL,
		indg2.opt=NULL,
		df0.opt=NULL,
		df1.opt=NULL,
		modtwo.opt=NULL)	
	}
}
