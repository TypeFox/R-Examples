#' @title Resuming the tree attributes: it gives general information about the tree
#' @details
#' Internal function. \code{mox.tree} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{mox.tree.pls}}. 
#' @return a data with general information about the tree and its nodes
#' @keywords internal
#' @export

mox.tree.pls	<-	function(tree,...)
{
	info.node		= list()
	type			= NULL
	terminal		= NULL
	perc			= NULL
	var			= NULL
	mox			= NULL
	if	(length(tree@nodes)>1)
	{
		for (n in tree@nodes)
		{
			if (n@id == 1)
			{
				length.root = length(n@elements)
			}
			if (length(n@childs) > 0)
			{
				info.node[[length(info.node)+1]] = data.frame(n@info@variable,n@id,n@childs,n@info@level)
			}
			if	(length(n@childs) == 0)
			{
				type		= "leaf"
				terminal	= "yes"
			}
			if	(n@father == 0)
			{
				type		= "root"
				terminal	= "no"
			}
			if	(n@father!=0 && length(n@childs) != 0)
			{
				type		= "node"
				terminal	= "no"	
			} 
			perc = round((length(n@elements)/length.root)*100,2)
			data = data.frame(n@id,n@father,showDeepth(n),type,terminal,length(n@elements),perc)		
			mox = rbind(mox,data)	
		}
	
		data.info.node = NULL
		
		for (i in 1:length(info.node)) {data.info.node = rbind(data.info.node,info.node[[i]])}	
	
		names(data.info.node)[2] = "n.father"
		names(data.info.node)[3] = "n.id"
	
		MOX =merge (mox, data.info.node,by="n.id",all.x=TRUE)[,-9]
	
		names(MOX) = c("Node","Parent","Depth","Type","Terminal","Size","Percent","Variable","Category")
	
		MOX	
	}
	else
	{
		MOX = NULL
	}
	MOX
}	
