#' @title Calculating the element of the root node
#' @details
#' Internal function. \code{root.tree} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element
#' @param \dots Further arguments passed on to \code{\link{root.tree}}. 
#' @return the elements of the root node
#' @keywords internal
#' @export

root.tree	<-	function(tree,...)
{
	root = NULL
	
	for (n in tree@nodes)
	{
		if (n@id == 1)
		{
			root=n@elements
		}
	}	
	root
}	
