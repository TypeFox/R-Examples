# wrapper over get.adjacency / get.edgelist


#'Convert igraph objects to adjacency or edge list matrices
#'
#'Get adjacency or edgelist representation of the network stored as
#'\code{igraph} object.
#'
#'If \code{matrix.type} is "edgelist" a two-column numeric edgelist matrix is
#'returned.  The value of \code{attrname} is ignored.
#'
#'If \code{matrix.type} is "adjacency" then a square adjacency matrix is
#'returned. If \code{attrname} is \code{NULL} (default) the matrix is binary.
#'Otherwise \code{attrname} can be a name of edge attribute of \code{x}. In
#'that case the cells of the results are the values of that attribute.
#'
#'Other arguments passed through \code{...} are passed to either
#'\code{\link[igraph]{get.adjacency}} or \code{\link[igraph]{get.edgelist}}
#'depending on the value of \code{matrix.type}.
#'
#'@param x object of class igraph, the network
#'@param matrix.type character, type of matrix to return, currently "adjacency"
#'or "edgelist" are supported
#'@param attrname character, name of the edge attribute to use to fill in the
#'cells of the adjacency matrix
#'@param sparse logical, whether to return a sparse matrix
#'@param \dots other arguments to/from other methods
#'@return Depending on the value of \code{matrix.type} either a square
#'adjacency matrix or a two-column numeric matrix representing the edgelist.
#'@seealso \code{\link[igraph]{get.adjacency}},
#'\code{\link[igraph]{get.edgelist}}
#'@examples
#'
#'data(exIgraph)
#'as.matrix(exIgraph, "adjacency")
#'as.matrix(exIgraph, "edgelist")
#'# use edge attribute "label" 
#'as.matrix(exIgraph, "adjacency", sparse=FALSE, attr="label")
#'
as.matrix.igraph <- function(x, matrix.type=c("adjacency", "edgelist"), attrname=NULL, sparse=FALSE, ...)
{
  mt <- match.arg(matrix.type)
  switch(mt,
      adjacency = igraph::get.adjacency(graph=x, attr=attrname, sparse=sparse, ...),
      edgelist = igraph::get.edgelist(graph=x, ...) )
}
