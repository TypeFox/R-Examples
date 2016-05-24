orderEdge <- function(adjMat) {
  # all credit goes to Markus Kalisch (kalisch@stat.math.ethz.ch) and Diego Colombo who
  # originally implemented this. We borrow this function and
  # slightly modified so as to meet our need.
  #
  # Kalisch, M., Machler, M., Colombo, D., Maathuis, M. H., &
  # Buehlmann, P. (2012). Causal inference using graphical models
  # with the R package pcalg. Journal of Statistical Software, 47(11), 1-26.

  stopifnot(ggm::isAcyclic(adjMat))
  ordered.nodes <- topoOrder(adjMat)
  edge.df <- makeEdgeDf(adjMat)
  eOrder <- 0
  while (any(unOrdered <- is.na(edge.df$order))) {
    counter <- 0
    y <- NA
    found <- FALSE
    while (!found) {
      counter <- counter + 1
      node <- ordered.nodes[counter]
      nbr.nodes <- which(adjMat[, node] == 1)
      if (length(nbr.nodes) > 0) {
        unlabeled <- rep.int(FALSE, length(nbr.nodes))
        for (i in seq_along(nbr.nodes)) {
          x <- nbr.nodes[i]
          unlabeled[i] <-
            length(intersect(which(edge.df$xmin == min(node, x)
                                   & edge.df$xmax == max(node, x)),
                             which(unOrdered))) > 0
        }
        if (any(unlabeled)) {
          nbr.unlab <- nbr.nodes[unlabeled]
          tmp <- ordered.nodes[ordered.nodes %in% nbr.unlab]
          y <- tmp[length(tmp)]
          edge.df$order[edge.df$xmin == min(node, y) &
                          edge.df$xmax == max(node, y)] <- eOrder
          eOrder <- eOrder + 1
          found <- TRUE
        }
      }
    }
  }
  return(edge.df)
}
