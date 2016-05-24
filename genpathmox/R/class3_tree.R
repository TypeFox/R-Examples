#' tree class
#'
#' tree is a S4 class that contains info on the binary segmentation tree
#'
#' @name tree-class
#' @rdname tree-class
						   					   
setClass("tree",representation(		
id			=	"numeric",
nodes	 	=	"list",
info_tree	=	"list"
))
