#' Function to generate a subgraph of a direct acyclic graph (DAG) induced by the input annotation data
#'
#' \code{dcDAGannotate} is supposed to produce a subgraph induced by the input annotation data, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph" or "Onto" object, a list of the vertices containing annotation data, and the mode defining the paths to the root of DAG. The induced subgraph contains vertices (with annotation data) and their ancestors along with the defined paths to the root of DAG. The annotations at these vertices (including their ancestors) are also updated according to the true-path rule: a domain annotated to a term should also be annotated by its all ancestor terms.
#'
#' @param g an object of class "igraph" or \code{\link{Onto}}
#' @param annotations an object of class \code{\link{Anno}}, that is, the vertices/nodes for which annotation data are provided
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced subgraph, an object of class "igraph" or "Onto" (the same as input). In addition to the original attributes to nodes and edges, the return subgraph is also appended by new node attributes: "annotations", which contains a list of domains either as original annotations or inherited annotations; "IC", which stands for information content defined as negative 10-based log-transformed frequency of domains annotated to that term.}
#' }
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @importFrom dnet dDAGinduce dDAGlevel visDAG dDAGtermSim dDAGgeneSim
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcEnrichment}}, \code{\link{dcDAGdomainSim}}, \code{\link{dcConverter}}
#' @include dcDAGannotate.r
#' @examples
#' \dontrun{
#' # 1) load onto.GOMF (as 'Onto' object)
#' g <- dcRDataLoader('onto.GOMF')
#'
#' # 2) load SCOP superfamilies annotated by GOMF (as 'Anno' object)
#' Anno <- dcRDataLoader('SCOP.sf2GOMF')
#'
#' # 3) prepare for annotation data
#' # randomly select 5 terms vertices (and their annotation data)
#' annotations <- Anno[,sample(1:dim(Anno)[2], 5)]
#'
#' # 4) obtain the induced subgraph according to the input annotation data
#' # 4a) based on all possible paths (i.e. the complete subgraph induced)
#' dcDAGannotate(g, annotations, path.mode="all_paths", verbose=TRUE)
#' # 4b) based on shortest paths (i.e. the most concise subgraph induced)
#' dag <- dcDAGannotate(g, annotations, path.mode="shortest_paths", verbose=TRUE)
#'
#' # 5) color-code nodes/terms according to the number of annotations
#' if(class(dag)=='Onto') dag <- dcConverter(dag, from='Onto', to='igraph')
#' data <- sapply(V(dag)$annotations, length)
#' names(data) <- V(dag)$name
#' dnet::visDAG(g=dag, data=data, node.info="both")
#' }

dcDAGannotate <- function (g, annotations, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), verbose=TRUE)
{
    
    path.mode <- match.arg(path.mode)
    
    if(class(g)=="Onto"){
        ig <- dcConverter(g, from='Onto', to='igraph', verbose=F)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'Onto' object.\n")
    }
    
    if (class(annotations) != "Anno"){
        stop("The function must apply to 'Anno' object.\n")
    }
    
    D <- annoData(annotations)
    originAnnos <- sapply(1:ncol(D), function(j){
        names(which(D[,j]!=0))
    })
    names(originAnnos) <- termNames(annotations)
    
    ## check nodes in annotations
    if (is.list(originAnnos)){
        originNodes <- names(originAnnos)
        
        ind <- match(originNodes, V(ig)$name)
        nodes_mapped <- originNodes[!is.na(ind)]
        if(length(nodes_mapped)==0){
            stop("The input annotations do not contain terms matched to the nodes/terms in the input graph.\n")
        }
    }
    
    ## generate a subgraph of a direct acyclic graph (DAG) induced by terms in input annotations
    dag <- dnet::dDAGinduce(ig, originNodes, path.mode=path.mode)
    allNodes <- V(dag)$name
    
    ## create a new (empty) hash environment
    ## node2domain.HoH: 1st key (node/term), 2nd key (domain), value (origin/inherit)
    node2domain.HoH <- new.env(hash=T, parent=emptyenv())
    
    ## assigin original annotations to "node2domain.HoH"
    lapply(allNodes, function(node){
        e <- new.env(hash=T, parent=emptyenv())
        if(node %in% originNodes){
            sapply(originAnnos[[node]], function(domain){
                assign(as.character(domain), "origin", envir=e)
            })
        }
        assign(node, e, envir=node2domain.HoH)
    })
  
    ## get the levels list
    level2node <- dnet::dDAGlevel(dag, level.mode="longest_path", return.mode="level2node")
    ## build a hash environment from the named list "level2node"
    ## level2node.Hash: key (level), value (a list of nodes/terms)
    level2node.Hash <- list2env(level2node)
    nLevels <- length(level2node)
    for(i in nLevels:1) {
        currNodes <- get(as.character(i), envir=level2node.Hash, mode='character')

        ## get the incoming neighbors (excluding self) that are reachable (i.e. nodes from i-1 level)
        adjNodesList <- lapply(currNodes, function(node){
            neighs.in <- igraph::neighborhood(dag, order=1, nodes=node, mode="in")
            setdiff(V(dag)[unlist(neighs.in)]$name, node)
        })
        names(adjNodesList) <- currNodes

        ## push the domains from level i to level i - 1
        lapply(currNodes, function(node){
            ## get the domains from this node
            domainsID <- ls(get(node, envir=node2domain.HoH, mode='environment'))

            ## assigin inherit annotations to "node2domain.HoH"
            lapply(adjNodesList[[node]], function(adjNode){
                adjEnv <- get(adjNode, envir=node2domain.HoH, mode='environment')
                sapply(domainsID, function(domainID){
                    assign(domainID, "inherit", envir=adjEnv)
                })
            })
        })
        
        if(verbose){
            message(sprintf("\tAt level %d, there are %d nodes, and %d incoming neighbors.", i, length(currNodes), length(unique(unlist(adjNodesList)))), appendLF=T)
        }
        
    }

    ## append 'annotations' attributes to the graph
    node2domains <- as.list(node2domain.HoH)[allNodes]
    domain_annotations <- sapply(node2domains, function(node){
        #unlist(as.list(node))
        names(unlist(as.list(node)))
    })
    V(dag)$annotations <- domain_annotations

    ## append 'IC' attributes to the graph
    counts <- sapply(domain_annotations, length)
    IC <- -1*log10(counts/max(counts))
    V(dag)$IC <- IC
    
    if(class(g)=="Onto"){
        dag <- dcConverter(dag, from='igraph', to='Onto', verbose=F)
    }
    
    return(dag)
}