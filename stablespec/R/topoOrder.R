topoOrder <- function(adjMat) {
  # all credit goes to Markus Kalisch (kalisch@stat.math.ethz.ch) and Diego Colombo who
  # originally implemented this. We borrow this function and
  # slightly modified so as to meet our need.
  #
  # Kalisch, M., Machler, M., Colombo, D., Maathuis, M. H., &
  # Buehlmann, P. (2012). Causal inference using graphical models
  # with the R package pcalg. Journal of Statistical Software, 47(11), 1-26.


  if (!ggm::isAcyclic(adjMat)) {
    stop("The graph is not acyclic!")
  }
  n <- nrow(adjMat)
  nod <- 1:n
  indeg <- rep(0, n)
  up <- !adjMat[lower.tri(adjMat)]
  if (all(up)) {
    return(nod)
  }
  zero.indeg <- c()
  for (i in nod) {
    indeg[i] <- sum(adjMat[, i])
    if (indeg[i] == 0)
      zero.indeg <- c(i, zero.indeg)
  }
  s <- 1
  ord <- rep(0, n)
  while (length(zero.indeg) > 0) {
    v <- zero.indeg[1]
    zero.indeg <- zero.indeg[-1]
    ord[s] <- v
    s <- s + 1
    cs <- nod[adjMat[v, ] == 1]
    if (length(cs) == 0) {
      next
    }
    for (j in 1:length(cs)) {
      k <- cs[j]
      indeg[k] <- indeg[k] - 1
      if (indeg[k] == 0) {
        zero.indeg <- c(k, zero.indeg)
      }
    }
  }
  return(ord)
}
