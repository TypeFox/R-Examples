#' Function to generate a subgraph induced by given vertices and their k nearest neighbors
#'
#' \code{dNetInduce} is supposed to produce a subgraph induced by given vertices and its k nearest neighbors. The input is a graph of "igraph" or "graphNET" object, a list of the vertices of the graph, and a k value for finding k nearest neighbors for these vertices. The output is a subgraph induced by given vertices plus their k neighbours. The resultant subgraph inherits the class from the input one. The induced subgraph contains exactly the vertices of interest, and all the edges between them. 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param nodes_query the vertices for which the calculation is performed
#' @param knn an integeter specifying how many k steps are used to find the nearest neighbours of the given vertices. By default, knn is set to zero; it means no neighbors will be considered. When knn is 1, the immediate neighbors of the given vertices will be also considered for inducing the subgraph. The same is true when knn is 2, etc
#' @param remove.loops logical to indicate whether the loop edges are to be removed. By default, it sets to false
#' @param largest.comp logical to indicate whether the largest component is only retained. By default, it sets to true for the largest component being left
#' @param min.comp.size an integer specifying the minimum size of component that will be retained. This parameter only works when setting the false to keep the largest component. By default, it sets to 1 meaning all nodes will be retained
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced subgraph, an object of class "igraph" or "graphNEL"}
#' }
#' @note The given vertices plus their k nearest neighbors will be used to induce the subgraph.
#' @export
#' @seealso \code{\link{dNetInduce}}
#' @include dNetInduce.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) select the first 10 vertices as the query nodes
#' nodes_query <- V(g)[1:10]
#'
#' # 3) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, nodes_query, knn=0)
#'
#' # 4) produce the induced subgraph based on the nodes in query ane their immediate neighbours
#' subg <- dNetInduce(g, nodes_query, knn=1)

dNetInduce <- function(g, nodes_query, knn=0, remove.loops=F, largest.comp=T, min.comp.size=1) 
{
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if(is.null(knn)){
        knn <- 0
    }
    knn <- as.integer(knn)
    if(knn < 0){
        knn <- 0
    }
    
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }
    ind <- match(nodes_query, V(ig)$name)
    nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
    
    nei <- unique(unlist(neighborhood(ig, nodes=nodes_mapped, order=knn)))
    subg <- igraph::induced.subgraph(ig, vids=nei)
    
    if(remove.loops){
        subg <- igraph::simplify(subg, remove.loops=T)
    }
    
    if(largest.comp==T){
        clust <- igraph::clusters(subg)
        cid <- which.max(clust$csize)
        subg <- igraph::induced.subgraph(subg, V(subg)[clust$membership==cid])
    }else{
    	min.comp.size <- as.integer(min.comp.size)
    	if(min.comp.size>=1 & min.comp.size<=igraph::vcount(subg)){
			clust <- igraph::clusters(subg)
			cid <- which(clust$csize >= min.comp.size)
			subg <- igraph::induced.subgraph(subg, V(subg)[clust$membership %in% cid])
    	}
    }
    
    if(class(g)=="graphNEL"){
        subg <- igraph.to.graphNEL(subg)
    }

    return(subg)
}