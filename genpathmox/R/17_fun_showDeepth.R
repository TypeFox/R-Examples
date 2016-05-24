#' @title Calculating Deepth stop criterion  
#' @details
#' Internal function. \code{showDeepth} is called by \code{pls.pathmox}.
#' @param node id that identify a specicif node
#' @return Deepth of the tree 
#' @keywords internal
#' @export

showDeepth=function(node)
{
	return (trunc(log2(node@id)))
}
