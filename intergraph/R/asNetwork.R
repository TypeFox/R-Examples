#' Convert objects to class "network"
#' 
#' Convert objects to class "network"
#' 
#' This is a generic function which dispatches on argument \code{x}.  It creates
#' objects of class "network" from other R objects.
#' 
#' The method for data frames is inspired by the similar function in package
#' \pkg{igraph}: \code{\link[igraph]{graph.data.frame}}.  It assumes that first
#' two columns of \code{x} constitute an edgelist.  The remaining columns are
#' interpreted as edge attributes. Optional argument \code{vertices} allows for
#' including vertex attributes.  The first column is assumed to vertex id, the
#' same that is used in the edge list. The remaining colums are interpreted as
#' vertex attributes.
#' 
#' The method for objects of class "igraph" takes the network of that class and
#' converts it to data frames using \code{\link{asDF}}. The network is recreated
#' in class "network" using \code{asNetwork.data.frame}. The function currently
#' does not support bipartite "igraph" networks.
#' 
#' @param x an R object to be coerced, see Details for the description of
#' available methods
#' @param amap data.frame with attribute copy/rename rules, see
#' \code{\link{attrmap}}
#' @param directed logical, whether the created network should be directed
#' @param vertices NULL or data frame, optional data frame containing vertex
#' attributes
#' @param \dots other arguments from/to other methods
#' @return Object of class "network".
#' @seealso \code{\link[igraph]{graph.data.frame}}
#' 
#' \code{\link{asIgraph}} for conversion in the other direction.
#'
#' @export
#'
#' @example examples/asNetwork.R
#'

asNetwork <- function(x, ...) UseMethod("asNetwork")

#' @method asNetwork data.frame
#' @export
#' @rdname asNetwork
asNetwork.data.frame <- function(x, directed=TRUE, vertices=NULL, ...)
{
  edb <- validateEL(x)
  # got vertex DB?
  if(!is.null(vertices))
  {
    vdb <- validateVDB(vertices)
    stopifnot(validNetDB(edb, vdb))
  }
  # number of vertices
  if(is.null(vertices)) nv <- length(unique(c(edb[,1], edb[,2])))
  else nv <- nrow(vertices)
  # create an empty network object
  rval <- network::network.initialize(nv, directed=directed, hyper=FALSE,
                                      multiple=any(duplicated(edb[,1:2])),
                                      loops=any(edb[,1] == edb[,2]))
  # add edges
  rval <- network::add.edges(rval, as.list(edb[,1]), as.list(edb[,2]))
  # add edge attribbutes
  if( ncol(edb) > 2)
    for(i in seq(3, ncol(edb)))
    {
      rval <- network::set.edge.attribute(rval, attrname=names(edb)[i], value=edb[,i])
    }
  # vertex attributes
  if( !is.null(vertices) && ncol(vertices) > 1 )
  {
    for( i in seq(2, ncol(vdb)) )
    {
      rval <- network::set.vertex.attribute(rval, attrname=names(vdb)[i],
                                            value=vdb[,i])
    }
  }
  rval
}


#' @method asNetwork igraph
#' @export
#' @rdname asNetwork
asNetwork.igraph <- function(x, amap=attrmap(), ...)
{
  object <- x
    na <- dumpAttr(object, "network")
    l <- asDF(object)

    ### prepare edge attributes
    eats <- attrmapmat("igraph", "network", "edge", db=amap)
    if( nrow(eats) > 0 )
    {
      # drop some
      todrop <- eats[ is.na(eats[,"toattr"]) , "fromattr" ]
      edges <- l$edges[ !( names(l$edges) %in% todrop ) ]
      # rename some
      names(edges) <- recode(names(edges), eats)
    } else
    {
      edges <-l$edges
    }

    ### prepare vertex attributes
    vats <- attrmapmat("igraph", "network", "vertex", db=amap)
    if( nrow(vats) > 0 )
    {
      # drop some
      todrop <- vats[ is.na(vats[,"toattr"]) , "fromattr" ]
      vertexes <- l$vertexes[ !( names(l$vertexes) %in% todrop )  ]
      # rename some
      names(vertexes) <- recode(names(vertexes), vats)
    } else
    {
      vertexes <- l$vertexes
    }

    ### make 'igraph' object
    rval <- asNetwork( edges,
        directed=igraph::is.directed(object),
        multiple = any(igraph::is.multiple(object)),
        loops = any(igraph::is.loop(object)),
        vertices=vertexes, ...)

    ### apply/rename/drop network attributes
    nats <- attrmapmat("igraph", "network", "network", db=amap)
    if( nrow(nats) > 0 )
    {
      todrop <- nats[ is.na(nats[,"toattr"]) , "fromattr" ]
      na <- na[ !( names(na) %in% todrop ) ]
      names(na) <- recode(names(na), nats)
    }
    if( length(na) > 0 )
    {
      for( naname in names(na) )
        network::set.network.attribute(rval, naname, na[[naname]])
    }
    if( is.function(network::get.network.attribute(rval, "layout")) )
      warning("network attribute 'layout' is a function, print the result might give errors")
    rval
}
