#' Function to find the tip node(s) of a direct acyclic graph (DAG)
#'
#' \code{dDAGtip} is supposed to find the tip node(s) of a direct acyclic graph (DAG; an ontology). It return the name (i.e Term ID) of the tip node(s). 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @return 
#' \itemize{
#'  \item{\code{tip}: the tip name (i.e. Term ID)}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dDAGtip}}
#' @include dDAGtip.r
#' @examples
#' \dontrun{
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) find tips
#' tips <- dDAGtip(g)
#' tips
#' }

dDAGtip <- function (g)
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
        ## get the outgoing neighbors (including self) that are reachable
        neighs.out <- igraph::neighborhood(ig, order=1, nodes=V(ig), mode="out")
        ## calculate the number of neighors (including self)
        num.neighs.out <- sapply(neighs.out, length)
        ## find the tips that does not have any outgoing neighbors (except self)
        tips <- V(ig)[num.neighs.out==1]$name
    }else{
        igr <- dDAGreverse(ig)
        tips <- dDAGroot(igr)
    }
    
    return(tips)
}