#' Function to reverse the edge direction of a direct acyclic graph (DAG)
#'
#' \code{dDAGreverse} is supposed to reverse the edge direction of a direct acyclic graph (DAG; an ontology). The return graph remains all attributes associated on nodes and edges. 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @return 
#' \itemize{
#'  \item{\code{gr}: a graph being reversed, an object of class "igraph" or "graphNEL"}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dDAGreverse}}
#' @include dDAGreverse.r
#' @examples
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) the graph with reverse edge direction
#' gr <- dDAGreverse(g)
#' gr

dDAGreverse <- function (g)
{

    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    ## get edge data frame with the first two columns ("from" and "to") and edge attributes (if any)
    e <- get.data.frame(ig, what="edges")
    ## swap "from" & "to"
    neworder <- 1:ncol(e)
    neworder[1:2] <- c(2,1)
    er <- e[neworder]
    colnames(er) <- colnames(e)
    ## create an igraph graph using ordered data
    gr <- graph.data.frame(d=er, vertices=get.data.frame(ig, what="vertices"), directed=TRUE)
    
    if(class(g)=="graphNEL"){
        gr <- igraph.to.graphNEL(gr)
    }

    return(gr)
}