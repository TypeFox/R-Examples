### This file contains functions for plot the unrooted tree.

plotnj <- function(unrooted.tree, X.class = NULL, type = "u", main = NULL,
    show.tip.label = FALSE, show.node.label = FALSE,
    edge.width = 1, edge.width.class = edge.width, ...){
  if(class(unrooted.tree) != "phylo"){
    stop("This is not a phylo tree.")
  }
  if(any(unrooted.tree$edge < 0)){
    stop("This is not a unrooted tree.")
  }
  if(is.null(X.class)){
    plot(unrooted.tree, type = type, show.tip.label = show.tip.label,
         main = main, ...)
  } else{
    unrooted.tree.edge.color <- rep("#000000", length(unrooted.tree$edge.length))
    unrooted.tree.edge.width <- rep(edge.width, length(unrooted.tree$edge.length))
    for(k in sort(unique(X.class))){
      unrooted.tree.edge.color[unrooted.tree$edge[, 2] %in% which(X.class == k)] <-
        .Color[(k - 1) %% length(.Color) + 1]
      unrooted.tree.edge.width[unrooted.tree$edge[, 2] %in% which(X.class == k)] <-
        edge.width.class
    }
    plot(unrooted.tree, type = type,
         edge.color = unrooted.tree.edge.color,
         edge.width = unrooted.tree.edge.width,
         show.tip.label = show.tip.label, main = main, ...)
  }
  if(show.node.label) nodelabels(frame = "n")
} # End of plotnj().
