#' Function to find heuristically maximum scoring subgraph
#'
#' \code{dNetFind} is supposed to find the maximum scoring subgraph from an input graph and scores imposed on its nodes. The input graph and the output subgraph are both of "igraph" or "graphNET" object. The input scores imposed on the nodes in the input graph can be divided into two parts: the positve nodes and the negative nodes. The searching for maximum scoring subgraph is deduced to find the connected subgraph containing the positive nodes as many as possible, but the negative nodes as few as possible. To this end, a heuristic search is used (see Note below).
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param scores a vector of scores. For each element, it must have the name that could be mapped onto the input graph. Also, the names in input "scores" should contain all those in the input graph "g", but the reverse is not necessary
#' @return 
#' a subgraph with a maximum score, an object of class "igraph" or "graphNEL"
#' @note The search procedure is heuristic to find the subgraph with the maximum score:
#' \itemize{
#' \item{i) transform the input graph into a new graph by collapsing connected positive nodes into a meta-node. As such, meta-nodes are isolated to each other but are linked via negative nodes (single-nodes). Clearly, meta-nodes have positive scores, and negative scores for the single-nodes.}
#' \item{ii) append the weight attribute to the edges in the transformed graph. There are two types of edges: 1) the single-single edge with two single-nodes as two ends, and 2) single-meta edge with a single-node as one end and a meta-node as the other end. The weight for a single-single edge is the absolute sum of the scores in its two-end single-nodes but normalised by their degrees. The weight for a single-meta edge is simply the absolute score in its single-node end normalised by the degree. As such, weights are all non-negative.}
#' \item{iii) find minimum spanning tree (MST) in the weighted transformed graph using Prim's greedy algorithm. A spanning tree of the weighted graph is a subgraph that is tree and connects all the node together. The MST is a spanning tree with the sum of its edge weights minimised amongst all possible spanning trees.}
#' \item{iv) find all shortest paths between any pair of meta-nodes in the MST. Within the weighted transformed graph in ii), a subgraph is induced containing nodes (only occuring in these shortest paths) and all edges between them.}
#' \item{v) within the induced subgraph, identify single-nodes that are direct neighbors of meta-nodes. For each of these single-nodes, also make sure it has the absolute scores no more than the sum of scores in its neighboring meta-nodes. These single-nodes meeting both criteria are called "linkers".}
#' \item{vi) still within the induced subgraph in v), find the linker graph that contains only linkers and edges between them. Similarly to iii), find MST of the linker graph, called 'linker MST'. Notably, this linker MST serves as the scaffold, which only contains linkers but has meta-nodes being direcly attached to.}
#' \item{vii) in linker MST plus its attached meta-nodes, find the optimal path that has the sum of scores of its nodes and attached meta-nodes maximised amongest all possible paths. Nodes along this optimal path plus their attached meta-nodes are called 'subgraph nodes'.}
#' \item{viii) finally, from the input graph extract a subgraph (called 'subgraph') that only contains subgraph nodes and edges betwen them. This subgraph is the maximum scoring subgraph containing the positive nodes as many as possible, but the negative nodes as few as possible.}
#' }
#' @export
#' @seealso \code{\link{dNetFind}}
#' @include dNetFind.r
#' @examples
#' # 1) generate an vector consisting of random values from beta distribution
#' x <- rbeta(1000, shape1=0.5, shape2=1)
#'
#' # 2) fit a p-value distribution under beta-uniform mixture model
#' fit <- dBUMfit(x, ntry=1, hist.bum=FALSE, contour.bum=FALSE)
#'
#' # 3) calculate the scores according to the fitted BUM and fdr=0.01
#' # using "pdf" method
#' scores <- dBUMscore(fit, method="pdf", fdr=0.05, scatter.bum=FALSE)
#' names(scores) <- as.character(1:length(scores))
#'
#' # 4) generate a random graph according to the ER model
#' g <- erdos.renyi.game(1000, 1/100)
#'
#' # 5) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 6) find the subgraph with the maximum score
#' subgraph <- dNetFind(subg, scores)

dNetFind <- function(g, scores)
{

    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }

    if(is.null(names(scores))) {
        stop("The function must require the names of the input scores.\n")
    }else if(any(is.na(names(scores)))){
        warning("Those scores with NA as names will be removed")
        scores <- scores[!is.na(names(scores))]
    }
    
    ## append the 'score' attribute to the vertices of the graph
    V(ig)$score <- scores[V(ig)$name]
    
    ## focus on those vertices with positive scores
    pos.nodes <- V(ig)[V(ig)$score>0]$name
    
    if(length(pos.nodes)==0){
        warning("No positive nodes")
        subgraph <- graph.empty(n=0, directed=F)
    }else if(length(pos.nodes)==1){
        subgraph <- dNetInduce(ig, pos.nodes, knn=0, remove.loops=T, largest.comp=T)
    }else{
        
        pos.subgraph <- dNetInduce(ig, pos.nodes, knn=0, remove.loops=T, largest.comp=F)
        
        ## Decompose a graph into components
        conn.comp.graph <- decompose.graph(pos.subgraph)
        ## Score each component
        score.comp <- unlist(lapply(lapply(conn.comp.graph,get.vertex.attribute,"score"), sum))
        ## Sort according to component scores
        ind_order <- order(score.comp, decreasing=T)
        conn.comp.graph <- conn.comp.graph[ind_order]
        score.comp <- score.comp[ind_order]
        ## For each component, append the score attribute with a single overall score to each graph
        for (i in 1:length(conn.comp.graph)) {
            conn.comp.graph[[i]]$score <- score.comp[i]
        }
        
        ######################################################
        ## create meta-nodes for connected component subgraphs
        v.id <- seq(1, vcount(ig))
        names(v.id) <- V(ig)$name
        
        ## Convert a graph to an edge list (without the name attribute)
        edgelist <- get.edgelist(ig, names=F)
        edgelist1 <- edgelist[,1]
        edgelist2 <- edgelist[,2]
        
        ## action for meta-nodes creation
        for (i in 1:length(conn.comp.graph)) {
            new.id <- vcount(ig) + i
            for (j in as.character(v.id[V(conn.comp.graph[[i]])$name])) {
                edgelist1[which(edgelist1==j)] <- new.id
                edgelist2[which(edgelist2==j)] <- new.id
            }
        }
        new.ids <- seq(vcount(ig)+1, vcount(ig)+length(conn.comp.graph))
        new.names <- paste("cluster", seq(1:length(conn.comp.graph)), sep="")
        names(new.ids) <- new.names
        
        v.id <- c(v.id, new.ids)
        v.name <- names(v.id)
        names(v.name) <- v.id
        new.edgelist <- cbind(v.name[as.character(edgelist1)], v.name[as.character(edgelist2)])
        
        ## create a graph containing meta-nodes
        mig <- graph.edgelist(new.edgelist, directed=F)
        
        ## make sure there are no-loops and multiple edges
        mig <- simplify(mig, remove.loops=T, remove.multiple=T)
        
        ######################################################
        ## append the score attribute to the nodes
        ### for nodes with negative scores
        node.score <- scores[V(mig)$name]
        names(node.score) <- V(mig)$name
        ### for meta-nodes
        node.score.cluster <- sapply(conn.comp.graph, get.graph.attribute, "score")
        names(node.score.cluster) <- new.names
        ind_cluster <- grep("cluster", names(node.score))
        node.score[ind_cluster] <- node.score.cluster[names(node.score[ind_cluster])]
        ### append
        V(mig)$score <- node.score
        
        ## append the score.degree attribute to the nodes
        ## score.degree is that zero for all positives, and absolute score normalised by degree for all negatives
        score.degree <- 1/(igraph::degree(mig) + 1)
        ### nodes with positive scores are replaced with zeros
        tmp_score <- V(mig)$score
        tmp_score[tmp_score>0] <- 0
        ### append
        V(mig)$score.degree <- score.degree * tmp_score
         
        ## append the weight attribute to the edges
        E(mig)$weight <- rep(0, length(E(mig)))
        tmp_n1 <- get.edgelist(mig,names=F)[,1]
        tmp_n2 <- get.edgelist(mig,names=F)[,2]
        E(mig)$weight <- -(V(mig)[tmp_n1]$score.degree + V(mig)[tmp_n2]$score.degree)
        
        ######################################################
        if(!is.connected(mig)){
            decomp.graphs <- decompose.graph(mig)
            sum.pos <- lapply(decomp.graphs, function(x) {
                sum(node.score[names(which(node.score[V(x)$name] > 0))])
            })
            mig <- decomp.graphs[[which.max(sum.pos)]]
            rm(decomp.graphs)
        }
        
        mst <- minimum.spanning.tree(mig, weights=E(mig)$weight)
        
        mst.cluster.id <- grep("cluster", V(mst)$name)
        names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
        ## sort meta-nodes (mst.cluster.id)
        tmp <- unlist(strsplit(names(mst.cluster.id),"cluster"))
        ttmp <- as.numeric(matrix(tmp, nrow=2)[2,])
        mst.cluster.id <- mst.cluster.id[order(ttmp)]
        
        all.ids <- c()
        # check if multiple clusters exist, 
        # if they do not exist, assign neg nodes to empty vector that is used later.
        if(length(mst.cluster.id)==1){
            neg.node.ids.2 = c()
        }else{
            ## pathes between each pair of meta-nodes, containing all meta-nodes and their linkers
            for (j in 1:(length(mst.cluster.id)-1)) {
                path <- get.all.shortest.paths(mst, from=mst.cluster.id[j], to=mst.cluster.id[(j+1):length(mst.cluster.id)])
                all.ids <- c(all.ids, unlist(path$res))
            }
            all.ids <- unique(all.ids)
            
            sub.mig <- dNetInduce(mig, V(mst)[all.ids]$name, knn=0, remove.loops=T, largest.comp=F)
            
            ## append the clusters attributes to the nodes (only those with negative scores)
            neg.node.ids <- which(V(sub.mig)$score<0)
            for (i in neg.node.ids) {
                tmp_nei <- neighbors(sub.mig, v=i)
                tmp_nei_meta <- grep("cluster",V(sub.mig)[tmp_nei]$name)
                V(sub.mig)[i]$clusters <- list(tmp_nei[tmp_nei_meta])
            }
               
            score.neg.nodes <- c()
            for (i in neg.node.ids) {
                if (!is.na(V(sub.mig)[i]$clusters[1])) {
                    borders <- c(i, V(sub.mig)[i]$clusters)
                    
                    ##########################
                    borders <- unlist(borders)
                    ##########################
                    
                    score.neg.nodes <- c(score.neg.nodes, sum(V(sub.mig)[borders]$score))
                }else{
                    score.neg.nodes <- c(score.neg.nodes, V(sub.mig)[i]$score)
                }
            }
            neg.node.ids.2 <- neg.node.ids[score.neg.nodes > 0]
        }
        
        if (length(neg.node.ids.2) == 0) {
            ## remain the meta-node with the highest score
            tmp <- unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)], "cluster"))
            ttmp <- as.numeric(matrix(tmp, nrow=2)[2,])
            tmp_nodes <- unlist(lapply(conn.comp.graph,get.vertex.attribute,"name")[ttmp])
            subgraph <- dNetInduce(ig, tmp_nodes, knn=0, remove.loops=F, largest.comp=T)
            
        }else{
            ## connected graph between borders (neg.node.ids.2)
            subg <- dNetInduce(sub.mig, V(sub.mig)[neg.node.ids.2]$name, knn=0, remove.loops=T, largest.comp=T)
            mst.subg <- minimum.spanning.tree(subg, E(subg)$weight)
        
            ## A function for getting path score
            getPathScore <- function(path, graph1, graph2){    
                ## score in graph1
                s1 <- V(graph1)[path]$score
                ## score in graph2
                tmp <- unique(unlist(V(graph1)[path]$clusters))
                s2 <- V(graph2)[tmp]$score
                ## return    
                sum(c(s1,s2))
            }
            
            ## find a best path in border graph (mst.subg) with the highest score
            max.score <- 0
            best.path <- c()
            for (i in 1:vcount(mst.subg)) {
                path <- get.all.shortest.paths(mst.subg, from=V(mst.subg)[i])
                path.score <- unlist(lapply(path$res, getPathScore, graph1=mst.subg, graph2=sub.mig))
                best.pos <- which.max(path.score)
                if (path.score[[best.pos]] > max.score) {
                    best.path <- path$res[[best.pos]]
                    max.score <- path.score[[best.pos]]
                }
            }
            
            if(length(best.path)!=1){
                cluster.list <- V(mst.subg)[best.path]$clusters
                names.list <- as.character(1:length(cluster.list))
                names(cluster.list) <- names.list
                names(best.path) <- names.list

                for (i in names.list) {
                    res <- lapply(cluster.list, intersect, cluster.list[[i]])
                    if(length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list)<as.numeric(i)))]), unlist(cluster.list[as.character(which(as.numeric(names.list)>as.numeric(i)))])))>0){
                        if (length(setdiff(res[[i]], unique(unlist(res[names(res)!=i]))))==0){
                            cluster.list <- cluster.list[names(cluster.list)!=i]
                            names.list <- names.list[names.list!=i]
                        }
                    }
                }
                best.path <- best.path[names.list]
            }
            
            ## extract those meta-nodes attached to best.path
            pos.cluster <- V(sub.mig)[unique(unlist(V(mst.subg)[best.path]$clusters))]$name
            tmp <- unlist(strsplit(pos.cluster, "cluster"))
            ttmp <- as.numeric(matrix(tmp, nrow=2)[2,])
            tmp_meta_nodes <- unlist(lapply(conn.comp.graph, get.vertex.attribute, "name")[ttmp])
            ## extract those border nodes along best.path
            tmp_border_nodes <- V(mst.subg)[best.path]$name
            
            tmp_nodes <- c(tmp_border_nodes, tmp_meta_nodes)
            subgraph <- dNetInduce(ig, tmp_nodes, knn=0, remove.loops=F, largest.comp=T)
        }
    }
    
    if(class(g)=="graphNEL"){
        subgraph <- igraph.to.graphNEL(subgraph)
    }
    
    return(subgraph)
}