#' node class
#'
#' node is a S4 class that contains info on each node of the binary segmentation tree
#'
#' @name node-class
#' @rdname node-class
				   					  					
setClass("node",representation(		
id			=	"numeric",
elements	=	"numeric",
father		=	"numeric",
childs		=	"numeric",
info		=	"info.pls"
))