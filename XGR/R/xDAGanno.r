#' Function to generate a subgraph of a direct acyclic graph (DAG) induced by the input annotation data
#'
#' \code{xDAGanno} is supposed to produce a subgraph induced by the input annotation data, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph", a list of the vertices containing annotation data, and the mode defining the paths to the root of DAG. The induced subgraph contains vertices (with annotation data) and their ancestors along with the defined paths to the root of DAG. The annotations at these vertices (including their ancestors) can also be updated according to the true-path rule: those annotated to a term should also be annotated by its all ancestor terms.
#'
#' @param g an object of class "igraph" to represent DAG
#' @param annotation the vertices/nodes for which annotation data are provided. It can be a sparse Matrix of class "dgCMatrix" (with variants/genes as rows and terms as columns), or a list of nodes/terms each containing annotation data, or an object of class 'GS' (basically a list for each node/term with annotation data)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to true
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced subgraph, an object of class "igraph". In addition to the original attributes to nodes and edges, the return subgraph is also appended by two node attributes: 1) "anno" containing a list of variants/genes either as original annotations (and inherited annotations; 2) "IC" standing for information content defined as negative 10-based log-transformed frequency of variants/genes annotated to that term.}
#' }
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xDAGanno.r
#' @examples
#' \dontrun{
#' # 1) SNP-based ontology
#' # 1a) ig.EF (an object of class "igraph" storing as a directed graph)
#' g <- xRDataLoader('ig.EF')
#'
#' # 1b) load GWAS SNPs annotated by EF (an object of class "dgCMatrix" storing a spare matrix)
#' anno <- xRDataLoader(RData='GWAS2EF')
#'
#' # 1c) prepare for annotation data
#' # randomly select 5 terms/vertices (and their annotation data)
#' annotation <- anno[, sample(1:dim(anno)[2],5)]
#' 
#' # 1d) obtain the induced subgraph according to the input annotation data
#' # based on shortest paths (i.e. the most concise subgraph induced)
#' dag <- xDAGanno(g, annotation, path.mode="shortest_paths", verbose=TRUE)
#'
#' # 1e) color-code nodes/terms according to the number of annotations
#' data <- sapply(V(dag)$anno, length)
#' names(data) <- V(dag)$name
#' dnet::visDAG(g=dag, data=data, node.info="both")
#' 
#' ####################
#' # Below is for those SNPs annotated by the term called 'ankylosing spondylitis'
#' # The steps 1a) and 1b) are the same as above
#' # 1c') prepare for annotation data
#' # select a term 'ankylosing spondylitis'
#' terms <- V(g)$term_id[grep('ankylosing spondylitis',V(g)$term_name, perl=TRUE)]
#' ind <- which(colnames(anno) %in% terms)
#' annotation <- lapply(ind, function(x){names(which(anno[,x]!=0))})
#' names(annotation) <- colnames(anno)[ind]
#' 
#' # 1d') obtain the induced subgraph according to the input annotation data
#' # based on all possible paths (i.e. the complete subgraph induced)
#' dag <- xDAGanno(g, annotation, path.mode="all_paths", verbose=TRUE)
#'
#' # 1e') color-code nodes/terms according to the number of annotations
#' data <- sapply(V(dag)$anno, length)
#' names(data) <- V(dag)$name
#' dnet::visDAG(g=dag, data=data, node.info="both")
#'
#' ###########################################################
#' # 2) Gene-based ontology
#' # 2a) ig.MP (an object of class "igraph" storing as a directed graph)
#' g <- xRDataLoader('ig.MP')
#'
#' # 2b) load human genes annotated by MP (an object of class "GS" containing the 'gs' component)
#' GS <- xRDataLoader(RData='org.Hs.egMP')
#' anno <- GS$gs # notes: This is a list
#'
#' # 2c) prepare for annotation data
#' # randomly select 5 terms/vertices (and their annotation data)
#' annotation <- anno[sample(1:length(anno),5)]
#' 
#' # 2d) obtain the induced subgraph according to the input annotation data
#' # based on shortest paths (i.e. the most concise subgraph induced) 
#' # but without applying true-path rule
#' dag <- xDAGanno(g, annotation, path.mode="shortest_paths", true.path.rule=TRUE, verbose=TRUE)
#'
#' # 2e) color-code nodes/terms according to the number of annotations
#' data <- sapply(V(dag)$anno, length)
#' names(data) <- V(dag)$name
#' dnet::visDAG(g=dag, data=data, node.info="both")
#' }

xDAGanno <- function (g, annotation, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=TRUE, verbose=TRUE)
{
    
    path.mode <- match.arg(path.mode)
    
    ig <- g
    if (class(ig) != "igraph"){
        stop("The function must apply to the 'igraph' object.\n")
    }
    
    if(class(annotation)=="GS"){
        originAnnos <- annotation$gs
    }else if(class(annotation)=="list"){
        originAnnos <- annotation
    }else if(class(annotation)=="dgCMatrix"){
		D <- annotation
		originAnnos <- sapply(1:ncol(D), function(j){
			names(which(D[,j]!=0))
		})
		names(originAnnos) <- colnames(annotation)
    }else{
    	stop("The input annotation must be either 'GS' or 'list' or 'dgCMatrix' object.\n")
    }
    
    ## check nodes in annotation
    if (is.list(originAnnos)){
        originNodes <- names(originAnnos)
        
        ind <- match(originNodes, V(ig)$name)
        nodes_mapped <- originNodes[!is.na(ind)]
        if(length(nodes_mapped)==0){
            stop("The input annotation data do not contain terms matched to the nodes/terms in the input graph.\n")
        }
    }
    
    ## generate a subgraph of a direct acyclic graph (DAG) induced by terms from input annotations
    dag <- dnet::dDAGinduce(ig, originNodes, path.mode=path.mode)
    allNodes <- V(dag)$name
    
	## create a new (empty) hash environment
	## node2domain.HoH: 1st key (node/term), 2nd key (domain), value (origin/inherit)
	node2domain.HoH <- new.env(hash=T, parent=emptyenv())
	
	## assigin original annotations to "node2domain.HoH"
	lapply(allNodes, function(node){
		e <- new.env(hash=T, parent=emptyenv())
		if(node %in% originNodes){
			sapply(originAnnos[[node]], function(x){
				assign(as.character(x), "origin", envir=e)
			})
		}
		assign(node, e, envir=node2domain.HoH)
	})

    ## whether true-path rule will be applied
    if(true.path.rule){
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
	}
	
	node2domains <- as.list(node2domain.HoH)[allNodes]
	domain_annotation <- lapply(node2domains, function(node){
		#unlist(as.list(node))
		names(unlist(as.list(node)))
	})
	
    ## append 'anno' attributes to the graph
    V(dag)$anno <- domain_annotation

    ## append 'IC' attributes to the graph
    counts <- sapply(V(dag)$anno, length)
    IC <- -1*log10(counts/max(counts))
    ### force those 'Inf' to be 'zero'
    if(1){
    	IC[is.infinite(IC)] <- 0
    }
    V(dag)$IC <- IC
    
    return(dag)
}
