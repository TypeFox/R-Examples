#' Function to find the root node of a direct acyclic graph (DAG)
#'
#' \code{dDAGroot} is supposed to find the root node of a direct acyclic graph (DAG; an ontology). It return the name (i.e Term ID) of the root node. 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @return 
#' \itemize{
#'  \item{\code{root}: the root name (i.e. Term ID)}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dDAGroot}}
#' @include dDAGroot.r
#' @examples
#' \dontrun{
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) find the root
#' root <- dDAGroot(g)
#' root
#' }

dDAGroot <- function (g)
{

    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    fast <- TRUE
    if(fast){
        ## get the incoming neighbors (including self) that are reachable
        neighs.in <- igraph::neighborhood(ig, order=1, nodes=V(ig), mode="in")
        ## calculate the number of neighors (including self)
        num.neighs.in <- sapply(neighs.in, length)
        ## find the root that does not have any incoming neighbors (except self)
        root <- V(ig)[num.neighs.in==1]$name
    }else{
        ## get edge data frame with the first two columns ("from" and "to") and edge attributes (if any)
        e <- get.data.frame(ig, what="edges")
        ## calculate how many times being a child for each node
        allnodes <- V(ig)$name
        times.child <- sapply(allnodes, function(child) sum(e[,2]==child))
        ## find the root that is never as a child
        root <- names(which(times.child==0))
    }
    
    return(root)
}