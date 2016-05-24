#' Dump network attributes to a list
#' 
#' Given a network return a list of all the attributes.
#' 
#' 
#' @aliases dumpAttr dumpAttr.network dumpAttr.igraph
#' @param x network object
#' @param type character, type of attributes to dump
#' @param \dots other arguments from/to other methods
#' @return If \code{type} is one of "network", "vertex" or "edge" then a list of
#' corresponding attributes.
#' 
#' If \code{type} is "all" then lists of lists of attributes.
#'
#' @export
#' @example examples/dumpAttr.R

dumpAttr <- function(x, ...) UseMethod("dumpAttr")

#' @method dumpAttr network
#' @export
#' @rdname dumpAttr
dumpAttr.network <- function(x, type=c("all", "network", "vertex", "edge"), ...)
{
	type <- match.arg(type)
	if(type == "all")
	{
		n <- c("network", "vertex", "edge")
		rval <- lapply( n, function(nam) dumpAttr(x, type=nam))
		names(rval) <- n
		return(rval)
	} else
	{
		type <- match.arg(type)
		nam <- switch(type,
		network = network::list.network.attributes(x),
		edge = network::list.edge.attributes(x),
		vertex = network::list.vertex.attributes(x) )
		rval <- switch( type,
		network = lapply( nam, function(a) network::get.network.attribute(x, a)),
		edge = lapply( nam, function(a) network::get.edge.attribute(x$mel, a)),
		vertex = lapply( nam, function(a) network::get.vertex.attribute(x, a)) )
		names(rval) <- nam
		return(rval)
	}
}


#' @method dumpAttr igraph
#' @export
#' @rdname dumpAttr
dumpAttr.igraph <- function(x, type=c("all", "network", "vertex", "edge"), ...)
{
	type <- match.arg(type)
	if(type == "all")
	{
		n <- c("network", "vertex", "edge")
		rval <- lapply( n, function(nam) dumpAttr(x, type=nam))
		names(rval) <- n
		return(rval)
	} else
	{
		nams <- switch( type,
			network = igraph::list.graph.attributes(x),
			edge = igraph::list.edge.attributes(x),
			vertex = igraph::list.vertex.attributes(x) )
		rval <- switch( type,
			network = lapply( nams, function(a) igraph::get.graph.attribute(x, a) ),
			edge = lapply( nams, function(a) igraph::get.edge.attribute(x, a) ),
			vertex = lapply( nams, function(a) igraph::get.vertex.attribute(x, a) ) )
		names(rval) <- nams
		return(rval)
	}
}
