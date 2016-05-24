#' Convert network to data frame(s)
#' 
#' Convert a network data object to, possibly two, data frames: a data frame
#' with an edge list with edge attributes (if any), and a data frame of vertexes
#' with vertex attributes (if any). This is a generic function, see below for
#' available methods.
#' 
#' Currently there are methods for \code{object} being in one of the following
#' classes: "network", "igraph".
#' 
#' The function first gets the graph edge list using the appropriate function
#' depending on the class of \code{object} (see below).  Edge attributes, if
#' any, are then extracted using \code{\link{dumpAttr}} and added to it.
#' 
#' The vertex data frame is constructed with a vertex id as a sequence of
#' integer numbers. Details are method-specific, see below.  Vertex attributes
#' are extracted with \code{\link{dumpAttr}} and added to this data frame.
#'
#' Method-specific notes:
#' 
#' @param object R object representing a network, see below for available
#' methods
#'
#' @param \dots other arguments passed to/from other methods
#'
#' @return List with two components:
#' \describe{
#' \item{\code{edges}}{containing an edge list data frame at first two columns
#' and edge attributes on further ones.}
#' \item{\code{vertexes}}{with vertex id in the first column, named \code{id}
#' and any vertex attributes in the other columns.}
#' }
#'
#' @export
#'
#' @example examples/asDF.R
#'
asDF <- function(object, ...) UseMethod("asDF")





#' @method asDF network
#' @export
#' @rdname asDF
#' @details
#' For objects of class "network". Objects of this class store the vertex ids
#' as integer numbers. There is also an attribute "vertex.names" which is
#' always created when using graph constructors provided in the package
#' \pkg{network}.  \code{asDF} adds "vertex.names" to the vertex data frame as
#' a normal attribute and does not use it as a vertex id in the edge list.
#' 
#' The edge list is created using \code{\link[network]{as.matrix.network}}
#' function and contains integer vertex ids.
asDF.network <- function(object, ...)
{
	# get edge list and substitute vertex names
    dfedge <- as.data.frame(network::as.matrix.network(object, "edgelist"),
        stringsAsFactors=FALSE)
	# add edge attributes, if any
    eattr <- dumpAttr(object, "edge")
    if( length(eattr) > 0 )
        dfedge <- cbind(dfedge, as.data.frame(eattr, stringsAsFactors=FALSE))
    # make vertex data frame
    dfvertex <- data.frame(intergraph_id=seq(1, network::network.size(object)))
    # add vertex attributes if any
    vattr <- dumpAttr(object, "vertex")
    if( length(vattr) > 0 )
        dfvertex <- cbind( dfvertex, as.data.frame(vattr, stringsAsFactors=FALSE))
    list(edges=dfedge, vertexes=dfvertex)
}





#' @method asDF igraph
#' @export
#' @rdname asDF
#' @details
#' Objects of class "igraph", as provided by the \pkg{igraph} package. Vertex
#' ids in these objects integers starting from 1 (in \pkg{igraph} version prior
#' to 0.6-0 vertex ids started from 0). However, it is also possible to provide
#' a vertex attribute "name". It is added to the vertex data frame as a normal
#' vertex attribute and is not used on the edge list data frame.
#' 
#' The edge list is created using \code{\link[igraph]{get.edgelist}} function
#' with argument \code{names} set to \code{FALSE} so that integer vertex ids
#' are used.
asDF.igraph <- function(object, ...)
{
    # get edgelist
    dfedge <- as.data.frame(igraph::get.edgelist(object, names=FALSE), 
        stringsAsFactors=FALSE)
	# add edge attributes, if any
    eattr <- dumpAttr(object, "edge")
    if( length(eattr) > 0 )
        dfedge <- cbind(dfedge, as.data.frame(eattr, stringsAsFactors=FALSE))
    # make vertex data frame
    dfvertex <- data.frame(intergraph_id=seq(1, igraph::vcount(object) ))
    # add vertex attributes, if any
    vattr <- dumpAttr(object, "vertex")
    if( length(vattr) > 0 )
        dfvertex <- cbind( dfvertex, as.data.frame(vattr, stringsAsFactors=FALSE))
    list(edges=dfedge, vertexes=dfvertex)
}
