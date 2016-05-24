#' Function to define/calculate the level of nodes in a direct acyclic graph (DAG)
#'
#' \code{dDAGlevel} is supposed to calculate the level of nodes, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph" or "graphNET" object, and the definition of the node level. The return can be the level for each node or the nodes for each level.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param level.mode the mode of how to define the level of nodes in DAG. It can be "longest_path" for defining the node level as the length of the longest path from the node to the root, and "shortest_paths" for defining the node level as the length of the shortest path from the node to the root
#' @param return.mode the mode of how to return the node level information. It can be "node2level" for returning a named vector (i.e. the level for each node), and "level2node" for returning a named list (i.e. nodes for each level)
#' @return 
#' When "return.mode" is "node2level", it returns a named vector: for each named node (i.e. Term ID), it stores its level
#' When "return.mode" is "level2node", it returns a named list: for each named level, it contains the names (i.e. Term ID) of nodes belonging to this level
#' @note 
#' The level for the root is 1. The level based on the longest path will ensure that nodes at the same level will never be reachable (i.e. in the same path), while the level based on the shortest path will not be necessary. The "longest path" based level can be useful in visiting nodes from the tipmost level to the root: 1) for the current node, all children have been visited; 2) nodes at the same level can be looked at independantly. The "shortest path" based level can be useful in deriving nodes according to their closeness to the root.
#' @export
#' @seealso \code{\link{dDAGroot}}, \code{\link{dDAGreverse}}
#' @include dDAGlevel.r
#' @examples
#' # 1) load HPPA as igraph object
#' ig.HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#'
#' # 2) randomly select vertices as the query nodes
#' nodes_query <- sample(V(g),5)$name
#'
#' # 3) obtain the complete subgraph induced
#' subg <- dDAGinduce(g, nodes_query)
#'
#' # 4) calculate the node levels
#' # 4a) definition based on the longest path
#' dDAGlevel(subg, level.mode="longest_path")
#' # 4b) definition based on the shortest path
#' dDAGlevel(subg, level.mode="shortest_path")
#' # 4c) definition based on the longest path, and return nodes for each level
#' dDAGlevel(subg, level.mode="longest_path", return.mode="level2node")

dDAGlevel <- function (g, level.mode=c("longest_path","shortest_path"), return.mode=c("node2level","level2node"))
{
    
    level.mode <- match.arg(level.mode)
    return.mode <- match.arg(return.mode)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    root <- dDAGroot(ig)
    if(length(root) >1){
        warning("The input DAG has multiple roots; recalculate the root after being reversed.\n")
        ig <- dDAGreverse(ig)
        root <- dDAGroot(ig)
    }
    if(is.null(root)){
        stop("The function must have the root; check the eligibility of the input DAG.\n")
    }else if(length(root) >1){
        stop("Even after automatic reversing, the input DAG still has multiple roots; check the eligibility of the input DAG.\n")
    }
    
    ## DAG being induced from nodes in query
    if(level.mode=="longest_path"){
        edgelist <- get.data.frame(ig, what="edges")
        
        ## create a new (empty) environment
        nodeLookDown <- new.env(hash=T, parent=emptyenv())
        
        ## initiate the level at the root
        parents <- root
        level <- 1
        
        ## iterate from the root to the tips
        while(length(parents)>0) {
            ## assign the current level
            sapply(parents, function(parent){
                assign(parent, level, envir=nodeLookDown)
            })
            ## get direct children
            children <- sapply(parents, function(parent){
                edgelist[edgelist[,1]==parent, 2]
            })
            ## move to the next level
            level <- level+1
            parents <- unique(unlist(children, use.names=F))
        }

        ## get node level
        node2level <- unlist(as.list(nodeLookDown))
        ## make sure it sorts according to V(ig)$name
        node2level <- node2level[V(ig)$name]
        
    }else if(level.mode=="shortest_path"){
        vpaths <- get.shortest.paths(ig, from=root, to=V(ig), output="vpath")
        if(length(vpaths)!=length(V(ig)$name)){
            vpaths <- vpaths$vpath
        }
        node2level <- sapply(1:length(vpaths), function(i) length(vpaths[[i]]))
        names(node2level) <- V(ig)$name
    }

    if(return.mode=="node2level"){
        return(node2level)
    }else if(return.mode=="level2node"){
        ## get a list of nodes for each level
        lvs <- sort(unique(node2level))
        level2node <- sapply(lvs, function(x) names(node2level[node2level==x]))
        names(level2node) <- lvs
        
        return(level2node)
    }
}