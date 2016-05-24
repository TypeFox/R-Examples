

#' @method as.igraph editmatrix
#' @param x An object of class \code{\link{editmatrix}}, \code{\link{editarray}} or \code{\link{editset}}
#' @param weighted see \code{\link[igraph]{graph.adjacency}}
#'
#' @export
#' @rdname adjacency
as.igraph.editmatrix <- function(x, nodetype=c("all", "rules","vars"), rules=editnames(x), vars=getVars(x), weighted=TRUE, ...){
    nodetype <- match.arg(nodetype)
    a <- adjacency(E=x, nodetype=nodetype, rules=rules, vars=vars, ...)
    g <- igraph::graph.adjacency(
          a, 
          weighted=weighted,
          mode = 'undirected'
        )
    # $type is handy for bipartite graph function in igraph...
    igraph::V(g)$type <- igraph::V(g)$vars <- attr(a, "vars")
    g
}

#' @method as.igraph editarray
#' @export
#' @rdname adjacency
as.igraph.editarray <- function(x, nodetype=c("all", "rules","vars"), rules=editnames(x), vars=getVars(x),weighted=TRUE, ...){
    nodetype <- match.arg(nodetype)
    a <- adjacency(E=x, nodetype=nodetype, rules=rules, vars=vars, ...)
    g <- igraph::graph.adjacency(
          a, 
          weighted=weighted,
          mode = 'undirected'
        )
    igraph::V(g)$type <- igraph::V(g)$vars <- attr(a, "vars")
    g
}


#' @method as.igraph editset
#' @export
#' @rdname adjacency
as.igraph.editset <- function(x, nodetype=c("all", "rules","vars"), rules=editnames(x), vars=getVars(x),weighted=TRUE, ...){
    nodetype <- match.arg(nodetype)
    a <- adjacency(E=x, nodetype=nodetype, rules=rules, vars=vars, ...)
    g <- igraph::graph.adjacency(
          a, 
          weighted=weighted,
          mode = 'undirected'
        )
    igraph::V(g)$type <- igraph::V(g)$vars <- attr(a, "vars")
    g
}
