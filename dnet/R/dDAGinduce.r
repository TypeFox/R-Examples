#' Function to generate a subgraph of a direct acyclic graph (DAG) induced by given vertices
#'
#' \code{dDAGinduce} is supposed to produce a subgraph induced by given vertices, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph" or "graphNET" object, a list of the vertices of the graph, and the mode defining the paths to the root of DAG. The resultant subgraph inherits the class from the input one. The induced subgraph contains exactly the vertices of interest and their defined paths to the root of DAG. 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param nodes_query the vertices for which the calculation is performed
#' @param path.mode the mode of paths induced by nodes in query. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced subgraph, an object of class "igraph" or "graphNEL"}
#' }
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{dDAGroot}}
#' @include dDAGinduce.r
#' @examples
#' \dontrun{
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) randomly select vertices as the query nodes
#' # the query nodes can be igraph vertex sequences
#' nodes_query <- sample(V(g),5)
#' # more commonly, the query nodes can be term id
#' nodes_query <- sample(V(g),5)$name
#'
#' # 3) obtain the induced subgraph
#' # 3a) based on all possible paths (i.e. the complete subgraph induced)
#' subg <- dDAGinduce(g, nodes_query, path.mode="all_paths")
#' # 3b) based on shortest paths (i.e. the most concise subgraph induced)
#' subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
#' }

dDAGinduce <- function (g, nodes_query, path.mode=c("all_paths","shortest_paths","all_shortest_paths"))
{
    
    path.mode <- match.arg(path.mode)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if(class(nodes_query)=="igraph.vs"){
        nodes_query <- nodes_query$name
    }
    
    ## check nodes in query
    ind <- match(nodes_query, V(ig)$name)
    nodes_mapped <- nodes_query[!is.na(ind)]
    if(length(nodes_mapped)==0){
        stop("Nodes in query cannot be found in the input graph.\n")
    }else{
        nodes_query <- V(ig)[nodes_mapped]
    }

    ## DAG being induced from nodes in query
    if(path.mode=="all_paths"){
        #edgelist <- get.data.frame(ig, what="edges")
        
        fast <- T
        if(fast){
            ## find all ancestors for any node
            neighs.in <- igraph::neighborhood(ig, order=vcount(ig), nodes=nodes_query$name, mode="in")
            nodeInduced <- V(ig)[unique(unlist(neighs.in))]$name
        }else{
            ## create a new (empty) hash environment: key (node), value (TRUE)
            node.Hash <- new.env(hash=T, parent=emptyenv())
            ## A function to iterate to the root, given a node
            buildInducedGraph <- function(node) {
                ## exists: true if and only if an object of the correct name and mode is found
                if (exists(node, envir=node.Hash, mode="logical", inherits=FALSE)){
                    ## for node already visited
                    return(1)
                }else{
                    ## assign the node (with the value 'TRUE') into node.Hash
                    assign(node, TRUE, envir=node.Hash)
                    ## get its direct parents
                
                    ## get the incoming neighbors (including self) that are reachable
                    neighs.in <- igraph::neighborhood(ig, order=1, nodes=node, mode="in")
                    adjNodes <- setdiff(V(ig)[unlist(neighs.in)]$name, node)
                    #adjNodes <- edgelist[edgelist[,2]==node, 1]
                    ## iterate until there are no direct parents
                    if (length(adjNodes)>0){
                        for (i in 1:length(adjNodes)){
                            buildInducedGraph(adjNodes[i])
                        }
                    }
                    ## finish lookup/loop
                    return(0)
                }
            }
            tmp <- lapply(nodes_query$name, buildInducedGraph)
            nodeInduced <- ls(node.Hash)
        }
            
    }else if(path.mode=="all_shortest_paths"){
        root <- dDAGroot(ig)
        aspaths <- get.all.shortest.paths(ig, from=root, to=nodes_query)
        nodeInduced <- unique(unlist(aspaths$res))
    }else if(path.mode=="shortest_paths"){
        root <- dDAGroot(ig)
        vpaths <- get.shortest.paths(ig, from=root, to=nodes_query, output="vpath")
        if(length(vpaths)!=length(nodes_query)){
            vpaths <- vpaths$vpath
        }
        
        nodeInduced <- unique(unlist(vpaths))
    }
    
    subg <- induced.subgraph(ig, vids=nodeInduced)
    
    if(class(g)=="graphNEL"){
        subg <- igraph.to.graphNEL(subg)
    }

    return(subg)
}