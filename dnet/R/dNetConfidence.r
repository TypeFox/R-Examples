#' Function to append the confidence information from the source graphs into the target graph
#'
#' \code{eConsensusGraph} is supposed to append the confidence information (extracted from a list of the source graphs) into the target graph. The confidence information is about how often a node (or an edge) in the target graph that can be found in the input source graphs. The target graph is an object of class "igraph" or "graphNEL", and the source graphs are a list of objects of class "igraph" or "graphNEL". It also returns an object of class "igraph" or "graphNEL"; specifically, the same as the input target graph but appended with the "nodeConfidence" attribute to the nodes and the "edgeConfidence" attribute to the edges.
#'
#' @param target the target graph, an object of class "igraph" or "graphNEL"
#' @param sources a list of the source graphs, each with an object of class "igraph" or "graphNEL". These source graphs will be used to calculate how often a node (or an edge) in the target graph that can be found with them. 
#' @param plot logical to indicate whether the returned graph (i.e. the target graph plus the confidence information on nodes and edges) should be plotted. If it sets true, the plot will display the returned graph with the size of nodes indicative of the node confidence (the frequency that a node appears in the source graphs), and with the width of edges indicative of the edge confidence (the frequency that an edge appears in the source graphs)
#' @return
#' an object of class "igraph" or "graphNEL", which is a target graph but appended with the "nodeConfidence" attribute to the nodes and the "edgeConfidence" attribute to the edges
#' @note None
#' @export
#' @import igraph supraHex
#' @seealso \code{\link{visNet}}
#' @include dNetConfidence.r
#' @examples
#' # 1) generate a target graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#' target <- dNetInduce(g, V(g), knn=0)
#'
#' # 2) generate a list source graphs according to the ER model
#' sources <- lapply(1:100, function(x) erdos.renyi.game(100*runif(1), 1/10))
#'
#' # 3) append the confidence information from the source graphs into the target graph
#' g <- dNetConfidence(target=target, sources=sources)
#'
#' # 4) visualise the confidence target graph
#' visNet(g, vertex.size=V(g)$nodeConfidence/10, edge.width=E(g)$edgeConfidence)

dNetConfidence <- function(target, sources, plot=F)
{
    
    ## check target
    if(class(target)=="graphNEL"){
        ig <- igraph.from.graphNEL(target)
    }else{
        ig <- target
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }
    
    ## check sources
    if(class(sources[[1]])=="graphNEL"){
        
        isources <- lapply(sources, function(x){
            igraph.from.graphNEL(x)
        })
        
    }else{
        isources <- sources
    }
    if (class(isources[[1]]) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(isources[[1]])$name)){
        for(i in 1:length(isources)){
            V(isources[[i]])$name <- as.character(V(isources[[i]]))
        }
    }
    
    #################################################
    ## union all nodes and edges in target and sources
    
    ## for sources
    nodes_in_sources <- lapply(isources, function(x) V(x)$name)
    edges_in_sources <- list()
    for(i in 1:length(isources)){
        is <- isources[[i]]
        tmp_edges <- get.edgelist(is)
        edges_in_sources[[i]] <- apply(tmp_edges, 1, function(x) paste(sort(x), collapse="|"))
    }
    
    ## for target
    nodes_in_target <- V(ig)$name
    edges_in_target <- apply(get.edgelist(ig), 1, function(x) paste(sort(x), collapse="|"))
    
    ## union
    nodes_union <- union(nodes_in_target, unlist(nodes_in_sources))
    edges_union <- union(edges_in_target, unlist(edges_in_sources))
    
    #################################################
    ## counts of nodes and edges in sources
    nodes_in_sources <- table(unlist(nodes_in_sources))
    edges_in_sources <- table(unlist(edges_in_sources))
    
    ## Those in sources falling into the union
    nnodes_union <- rep(0, length(nodes_union))
    names(nnodes_union) <- nodes_union
    nnodes_union[names(nodes_in_sources)] <- nodes_in_sources
    
    nedges_union <- rep(0, length(edges_union))
    names(nedges_union) <- edges_union
    nedges_union[names(edges_in_sources)] <- edges_in_sources
    
    ## how many nodes and edges in target are found in union
    nnodes_found <- nnodes_union[nodes_in_target]
    nedges_found <- nedges_union[edges_in_target]
    
    #########################################
    ## calculate the frequency in the measure of the confidence
    node.confidence <- nnodes_found/length(isources)
    edge.confidence <- nedges_found/length(isources)
    
    ## cope the input graph g and append the nodeConfidence attribute to the nodes and the edgeConfidence attribute to the edges
    icg <- ig
    V(icg)$nodeConfidence <- node.confidence * 100
    E(icg)$edgeConfidence <- edge.confidence * 100
    
    if(plot){
        if(0){
            rescale <- function(x){
                (x-min(x))/(max(x)-min(x))
            }
            tmp_size <- 10*rescale(V(icg)$nodeConfidence + 1)
            tmp_width <- 10*rescale(E(icg)$edgeConfidence + 1)
        }else{
            if(max(V(icg)$nodeConfidence)>15){
                tmp_size <- (1+V(icg)$nodeConfidence)/100
            }else{
                tmp_size <- (1+V(icg)$nodeConfidence)/10
            }
            
            if(max(E(icg)$edgeConfidence)>15){
                tmp_width <- (1+E(icg)$edgeConfidence)/10
            }else{
                tmp_width <- (1+E(icg)$edgeConfidence)/10
            }
        }
        
        visNet(icg, vertex.size=(1+tmp_size), edge.width=(1+tmp_width))
    }
    
    if(class(target)=="graphNEL"){
        cg <- igraph.to.graphNEL(icg)
    }else{
        cg <- icg
    }
    
    return(cg)
    
}