#' @title Calculating the element of the terminal nodes
#' @details
#' Internal function. \code{terminal.tree} is called by \code{pls.pathmox}.
#' @param tree class containing the tree element.
#' @param \dots Further arguments passed on to \code{\link{terminal.tree}}. 
#' @return the elements of the terminal node
#' @keywords internal
#' @export

terminal.tree	<- function(tree,...)
{
	terminal 	= list()
	id 			= list()

	if (length(tree@nodes) > 1)
	{
		for (n in tree@nodes)
		{
			if (n@id == 1)
			{
				terminal[[length(terminal)+1]] = n@elements
				id[[length(id)+1]] = "Root"
			}
			if (length(n@childs) == 0)
			{
				terminal[[length(terminal)+1]] = n@elements
				id[[length(id)+1]] = n@id
			}
		}
		for (i in 1:length(terminal)){names(terminal) = paste("node",id)}
	
		terminal
	}
	else
	{
		terminal = NULL
	}
	terminal
}	
