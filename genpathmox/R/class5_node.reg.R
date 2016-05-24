#' node.reg class
#'
#' info.pls is a S4 class that contains element of the node class
#'
#' @name node.reg_class
#' @rdname node.reg_class

setClass("node.reg",representation(
id				= "numeric",
elements		= "numeric",
father 	    	= "numeric",
childs			= "numeric",
info 			= "info.reg"
))
