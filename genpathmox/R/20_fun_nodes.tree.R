#' @title Calculating the element of all nodes
#' @details
#' Internal function. \code{nodes.tree} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element
#' @param \dots Further arguments passed on to \code{\link{nodes.tree}}. 
#' @return the elements of all node
#' @keywords internal
#' @export

nodes.tree	<-	function(tree,...)
{
	nodes 	= list()
	id 		= list()
	
	if (length(tree@nodes) > 1)
	{
		for (n in tree@nodes)
		{
			nodes[[length(nodes)+1]] = n@elements
			id[[length(id)+1]] = n@id
		}
		for (i in 1:length(nodes))	{names(nodes) = paste("node",id)}
	
		nodes
	}
	else
	{
		nodes = NULL
	}
	nodes
}	
