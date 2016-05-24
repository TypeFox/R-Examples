#' @title Listing the the f-coefficients test results for each obtained partition 
#' @details
#' Internal function. \code{fcoef.tree.pls} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{fcoef.tree.pls}}. 
#' @return the f-coefficients test results for each obtained partition
#' @keywords internal
#' @export
	
fcoef.tree.pls	<-	function(tree,...)
{
	fc		= list()
	id		= list()
	fctable	= NULL
	if (length(tree@nodes) > 1)
	{
		for (n in tree@nodes)
		{
			if (length(n@childs) > 0)
			{
				id[[length(id)+1]] = n@id
				fctable	=	data.frame(as.matrix(n@info@fcstatistic),as.matrix(n@info@fpvalc))
				colnames(fctable)	= c("fc.statistic","fc.pvalue")
				fc[[length(fc)+1]]	= fctable
			}
		}
		
		names(fc) = paste("node",id,sep="")
		Fc.r = fc
	}
	else
	{
		Fc.r=list(fc=NULL,Signif=NULL)	
	}
	Fc.r
}
