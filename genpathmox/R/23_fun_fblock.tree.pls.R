#' @title Listing the the f-block test results for each obtained partition 
#' @details
#' Internal function. \code{fblock.tree.pls} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{fblock.tree.pls}}. 
#' @return the f-block test results for each obtained partition
#' @keywords internal
#' @export

fblock.tree.pls	<-	function(tree,...)
{
	fblock		= list()
	pfblock		= list()
	id			= list()
	fbtable 		= NULL
	pfbtable 		= NULL
	
	if (length(tree@nodes) > 1)
	{
		for (n in tree@nodes)
		{
			if (length(n@childs) > 0 && length(n@info@fbstatistic) > 0)
			{
				fblock[[length(fblock)+1]] = n@info@fbstatistic
				pfblock[[length(pfblock)+1]] = n@info@fpvalb
				id[[length(id)+1]] = n@id
			}
		}
		for(i in 1:length(fblock)) 
		{
			fbtable = rbind(fbtable,fblock[[i]])
			pfbtable = rbind(pfbtable,pfblock[[i]])
		}
	
		rownames(fbtable)  = paste("node",unlist(id))
		rownames(pfbtable) = paste("node",unlist(id))
		
		Fb.r = list(fbtable=fbtable,pfbtable=pfbtable)
	}
	else
	{
		Fb.r = list(fbtable=NULL,pfbtable=NULL)
	}
	Fb.r
}	
