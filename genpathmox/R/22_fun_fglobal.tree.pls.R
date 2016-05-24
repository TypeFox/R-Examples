#' @title Listing the the f-global test results for each obtained partition 
#' @details
#' Internal function. \code{fglobal.tree} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{fglobal.tree.pls}}. 
#' @return the f-global test results for each obtained partition
#' @keywords internal
#' @export

fglobal.tree.pls	<-	function(tree,...)
{
	fglobal = list()
	fgtable = NULL
	
	if (length(tree@nodes) > 1)
	{
		for (n in tree@nodes)
		{
			if (length(n@childs) > 0)
			{
				fglobal[[length(fglobal)+1]] = data.frame(n@id,n@info@fgstatistic,n@info@fpvalg,n@info@variable,t(n@info@level))
			}
		}
	
		for (i in 1:length(fglobal)) {fgtable = rbind(fgtable,fglobal[[i]])}
	
		colnames(fgtable)	=	c("node","fg.statistic","fg.pvalue","variable","g1.mod","g2.mod")
	
		Fg.r = fgtable
	}
	else
	{
		Fg.r = NULL
	}
	Fg.r
}	
