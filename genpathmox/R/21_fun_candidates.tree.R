#' @title Listing the candidates partions for each intermedied nodes 
#' @details
#' Internal function. \code{nodes.tree} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{candidates.tree}}. 
#' @return the candidates partions for intermediated nodes 
#' @keywords internal
#' @export

candidates.tree	<-	function(tree,...)
{
	candidates = list()
	id = list()

	if (length(tree@nodes) > 1)
	{
		for (n in tree@nodes)
		{
			if (length(n@childs)>0)
			{
			candidates[[length(candidates)+1]] = n@info@candidates
			id[[length(id)+1]] = n@id
			}
		}
		for (i in 1:length(candidates))	{names(candidates) = paste("node",id)}
	
		candidates
	}
	else
	{
		candidates = NULL
	}
	candidates
}	
