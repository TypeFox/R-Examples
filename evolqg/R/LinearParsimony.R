#'Estimation of Linear Parcimony ancestral stated for multivariate continous traits
#'
#'Wagner optimization
#'
#'@param tree phylogenetic tree
#'@param tip_means list of tip nodes data. Expects one vector for each terminal.
#'@param ... aditional arguments for phytools::fastAnc
#'@return list with calculated ancestral states, using labels or numbers from tree
#'@importFrom ape reorder.phylo 
#'@importFrom stats median
#'@import plyr
#'@examples
#'
#'data(dentus)
#'data(dentus.tree)
#'mean.list <- dlply(dentus, .(species), numcolwise(mean))
#'names(mean.list) <- dentus.tree$tip.label
#'evolqg:::LinearParsimony(dentus.tree, mean.list)
LinearParsimony <- function(tree, tip_means, ...){
    if(is.null(tree$node.label)){
        node.names <- tree$tip.label
    } else{
        node.names <- c(tree$tip.label, tree$node.label)
    }
    ancestral.stats <- list()
    for(tip in 1:length(tip_means)){
        ancestral.stats[[names(tip_means)[tip]]] = as.numeric(tip_means[[tip]])
    }
    if(!all(tree$tip.label %in% names(tip_means))) stop("All tip labels must be in stat list.")

    #Down pass for initial medians
    node.order <- unique(reorder(tree, "postorder")$edge[,1])
    for (node in node.order){
        if(is.na(node.names[node])) node.names[node] <- as.character(node)
        descendants.list <- node.names[tree$edge[which(tree$edge[,1]==node),2]]
        ancestral.stats[[node.names[node]]] <- spatialMedian(ancestral.stats[descendants.list])
    }
    total_dist = totalDist(tree$edge, ancestral.stats, node.names)

    #First up pass using descendants and ancestral median
    for (node in rev(node.order)){
        descendants.list <- c(node.names[tree$edge[which(tree$edge[,1]==node), 2]],
                              node.names[tree$edge[which(tree$edge[,2]==node), 1]])
        ancestral.stats[[node.names[node]]] <- spatialMedian(ancestral.stats[descendants.list])
    }
    total_dist = totalDist(tree$edge, ancestral.stats, node.names)

    #Up pass using descendants and ancestral median
    for (node in (node.order)){
        descendants.list <- c(node.names[tree$edge[which(tree$edge[,1]==node), 2]],
                              node.names[tree$edge[which(tree$edge[,2]==node), 1]])
        ancestral.stats[[node.names[node]]] <- spatialMedian(ancestral.stats[descendants.list])
    }
    total_dist = totalDist(tree$edge, ancestral.stats, node.names)
    return(ancestral.stats)
}

#'@importFrom depth med
spatialMedian <- function(x) {
    if(length(x[[1]]) == 1)
        median(unlist(x))
    else
        suppressWarnings(med(ldply(x, identity, .id = NULL), method = "Spatial"))$median
}
euclidianDist <- function(x, y) sqrt(sum((x - y)**2))
totalDist <- function(edges, ancestral.stats, node.names) sum(aaply(edges, 1, function(edge) euclidianDist(ancestral.stats[[node.names[[edge[1]]]]], ancestral.stats[[node.names[[edge[2]]]]])))
