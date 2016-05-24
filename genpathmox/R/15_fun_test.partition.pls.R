#' @title Calculating the comparison tests  
#' @details
#' Internal function. \code{test.partition.pls} is called by \code{pls.pathmox}.
#' @param x matrix or data frame containing the data.
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
#' @param modtwo vector indicating the binary partition
#' @param signific string indicating a stop condition.
#' @param method string indicating the method: LM or LAD.
#' @param \dots Further arguments passed on to \code{\link{test.partition.pls}}. 
#' @return list containing matrices needed for the comparison test
#' @keywords internal
#' @export

test.partition.pls <- function(x,inner,outer,mode,scheme,scaling,scaled,modtwo,signif,method,...) 
{
	d.info	=	F.data.pls(x,inner,outer,mode,scheme,scaling,scaled,modtwo)
	
	FG		=	Fg.test.pls(d.info$Y0,d.info$X0,d.info$Y1,d.info$X1,method)

	if(FG$pvg > signif)
	{
		list(Fg=FG$Fg ,pvg=FG$pvg,Fb=NULL ,pvb=NULL,Fc=list(),pvc=list())
	}
	else
	{
		FB		=	Fb.test.pls(d.info$Y1,d.info$X1,d.info$info.block,method)
		
		if(length(FB$pvb[FB$pvb<=signif] != 0))
		{
			FC		=	Fc.test.pls(d.info$Y1,d.info$X1,d.info$path.name,d.info$info.block,method)

			list(Fg=FG$Fg ,pvg=FG$pvg,Fb=FB$Fb ,pvb=FB$pvb,Fc=FC$Fc,pvc=FC$pvc) 	
		}
		else
		{
			list(Fg=FG$Fg ,pvg=FG$pvg,Fb=FB$Fb ,pvb=FB$pvb,Fc=list(),pvc=list())
		}
	}
}

