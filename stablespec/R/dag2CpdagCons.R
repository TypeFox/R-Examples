dag2CpdagCons <- function(dag, consMatrix) {
  # all credit goes to Markus Kalisch (kalisch@stat.math.ethz.ch) and Diego Colombo who
  # originally implemented this. We borrow this function and
  # slightly modified so as to meet our need.
  #
  # Kalisch, M., Machler, M., Colombo, D., Maathuis, M. H., &
  # Buehlmann, P. (2012). Causal inference using graphical models
  # with the R package pcalg. Journal of Statistical Software, 47(11), 1-26.

  p <- graph::numNodes(dag)
  if (graph::numEdges(dag) == 0) {
    cpdag.res <- dag
  } else {
    dag <- as(dag, "matrix")
    dag[dag != 0] <- 1
    e.df <- toLabelEdge(dag, consMatrix)
    cpdag <- matrix(0, p, p)
    for (i in 1:dim(e.df)[1]) {
      if (e.df$label[i]) {
        cpdag[e.df$tail[i], e.df$head[i]] <- 1
      } else {
        cpdag[e.df$tail[i], e.df$head[i]] <- cpdag[e.df$head[i], e.df$tail[i]] <- 1
      }
    }
    rownames(cpdag) <- colnames(cpdag) <- as.character(seq(1, p))
    cpdag.res <- as(cpdag, "graphNEL")
  }
  return(cpdag.res)
}
