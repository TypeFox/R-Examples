#' @title Listing the the f-coefficients test results for regression model for each obtained partition 
#' @details
#' Internal function. \code{fcoef.tree.reg} is called by \code{reg.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{fcoef.tree.reg}}.
#' @return the f-coefficients test results for each obtained partition
#' @keywords internal
#' @export

fcoef.tree.reg	<-	function(tree,...)
{
	fc=list()
	id=list()
	fctable = NULL
	if(length(tree@nodes)>1)
	{
		for (n in tree@nodes)
		{
			if (length(n@childs)>0)
			{
				id[[length(id)+1]]=n@id
				sign.fcpval=sign(n@info@fpvalc)
				fctable=data.frame(as.matrix(n@info@fcstatistic),as.matrix(n@info@fpvalc))
				colnames(fctable)=c("fc.statistic","fc.pvalue")
				fc[[length(fc)+1]]=fctable
			}
		}
		
		names(fc)=paste("node",id,sep="")
		Fc.r=fc
	}
	else
	{
		Fc.r=list(fc=NULL,Signif=NULL)	
	}
	Fc.r
}

