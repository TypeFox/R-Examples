#' @title Listing the the f-global test results for regression model for each obtained partition 
#' @details
#' Internal function. \code{fglobal.tree.reg} is called by \code{reg.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{fglobal.tree}}.
#' @return the f-global test results for each obtained partition
#' @keywords internal
#' @export

fglobal.tree	<-	function(tree,...)
{
	fglobal=list()
	fgtable = NULL
	
	if(length(tree@nodes)>1)
	{
		for (n in tree@nodes)
		{
			if (length(n@childs)>0)
			{
				sign.fgpval=sign(n@info@fpvalg)
				fglobal[[length(fglobal)+1]]=data.frame(n@id,n@info@fgstatistic,n@info@fpvalg,n@info@variable,t(n@info@modalidad))
			}
		}
		for(i in 1:length(fglobal)) {fgtable=rbind(fgtable,fglobal[[i]])}
		
		colnames(fgtable)	=	c("node","fg.statistic","fg.pvalue","variable","mod.g1","mod.g2")
		
		Fg.r=fgtable	
	}
	else
	{
		Fg.r=NULL
	}
	Fg.r
}	
