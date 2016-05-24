#' @title printing the tree structure
#' @details
#' Internal function. \code{printTree} is called by \code{reg.pathmox}.
#' @param tree that identify the tree object
#' @return the tree structure 
#' @keywords internal
#' @export

printTree	<-	function(tree)
{
	for (n in tree@nodes){
		print (n)
	}
}
