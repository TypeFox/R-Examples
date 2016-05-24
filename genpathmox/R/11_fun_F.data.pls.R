#' @title Defining the matrices needed for the comparison test  
#' @details
#' Internal function. \code{F.data.pls} is called by \code{test.partition.pls}.
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
#' @param modtwo vector indicating the binary partition
#' @param \dots Further arguments passed on to \code{\link{F.data.pls}}. 
#' @return list containing matrices needed for the comparison test
#' @keywords internal
#' @export

F.data.pls	<-	function(x,inner,outer,mode,scheme,scaling,scaled,modtwo,...)
{

	pls	=	plspm(x,inner,outer,mode,scaling,scheme,scaled=scaled)

	
	pls 		=	plspm(x,inner,outer,mode,scaling,scheme,scaled=scaled)
	
	LV = pls$scores
	
	path.name		=	element(pls$path_coefs)
	
	global		=	build.block(inner=inner,latent=LV)

	Y0			=	global$resp
	X0			=	global$x.block
	info.block	=	global$constructs.label

	g1.latent = subset(LV,modtwo==1)                    
	g2.latent = subset(LV,modtwo==2)
	
	g1=build.block(inner,latent=g1.latent)
	g2=build.block(inner,latent=g2.latent)
		
	g1.resp	=	g1$resp
	g1.pred	=	g1$x.block

	g2.resp	=	g2$resp
	g2.pred	=	g2$x.block

	Y.alt	=	rbind(g1.resp,g2.resp)
	X.alt	=	blockdiag(g1.pred,g2.pred)
	
	n.col	=	ncol(X.alt)/2
	
	colnames(X.alt)[1:n.col]	=	paste("g1 -",colnames(X0),sep=" ")
	colnames(X.alt)[(n.col+1):(2*n.col)]	=	paste("g2 -",colnames(X0),sep=" ")
	
	list(Y0=Y0,
		X0=X0,
		Y1=Y.alt,
		X1=X.alt,
		path.name=path.name,
		info.block=info.block)	
}

