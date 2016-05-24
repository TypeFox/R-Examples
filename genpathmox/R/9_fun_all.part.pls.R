#' @title Defining the candidates to the optimum partition for each of segmentation variables  
#' @details
#' Internal function. \code{all.part.pls} is called by \code{partopt.pls}.
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
#' @param \dots Further arguments passed on to \code{\link{all.part.pls}}. 
#' @return list containing information of the candidates to the optimum partition for each of segmentation variables  
#' @keywords internal
#' @export

all.part.pls <- function(x,y,inner,outer,mode,scheme,scaling,scaled,method,n.node,...) 
{
	part			=	partition(y)
	                                 
	p.bin		=	list()                                     
	level		=	list()                                    
	modtwo		=	list() 
	                                   
	length(p.bin)	=	length(level)	=	length(modtwo)	=	ncol(y) 
	
	Ftest      	= NULL
	df0        	= NULL
	df1        	= NULL
	pvl        	= NULL 
	g1.ind    	= NULL
	g2.ind		= NULL
	
	for (j in 1:ncol(y))
	{                              
		if (is.null(part$split[[j]])) next
		
		jpart	=	splitopt.pls(x,inner,outer,mode,scheme,scaling,scaled,splits=part$split[[j]],fact=y[,j],method,n.node) 
		
		if (is.null(jpart$pval)) next
		
		p.bin[[j]]			=	jpart$partition                   
		level[[j]]			=	unlist(jpart$mod)
		pvl[j]				=	jpart$pval
		Ftest[j]				=	jpart$F
		g1.ind[j]				=	jpart$g1.ind
		g2.ind[j]				=	jpart$g2.ind
		modtwo[[j]]			=	jpart$modtwo
		df0[j]				=	jpart$df.num
		df1[j]				=	jpart$df.den
	}
	variable					=	names(y)
	
	list(p.bin=p.bin,
		variable=variable,
		level=level,
		Ftest=Ftest,
		pvl=pvl,
		g1.ind=g1.ind,
		g2.ind=g2.ind,
		df0=df0,
		df1=df1,
		modtwo=modtwo)	
}

