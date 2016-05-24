makeEdgeDf <- function(amat) {
  # all credit goes to Markus Kalisch (kalisch@stat.math.ethz.ch) and Diego Colombo who
  # originally implemented this. We borrow this function and
  # slightly modified so as to meet our need.
  #
  # Kalisch, M., Machler, M., Colombo, D., Maathuis, M. H., &
  # Buehlmann, P. (2012). Causal inference using graphical models
  # with the R package pcalg. Journal of Statistical Software, 47(11), 1-26.


  stopifnot(sum(amat) > 0)
  e <- which(amat == 1, arr.ind = TRUE)
  e.dup <- duplicated(t(apply(e, 1, sort)))
  nmb.edges <- sum(!e.dup)
  res <- data.frame(xmin = rep(NA, nmb.edges), xmax = rep(NA, nmb.edges),
                    tail = rep(NA, nmb.edges), head = rep(NA, nmb.edges),
                    order = rep(NA, nmb.edges), type = rep(1, nmb.edges))

  pure.edges <- e[!e.dup, ]
  if (length(pure.edges) == 2) {
    dim(pure.edges) <- c(1, 2)
  }
  for (i in 1:dim(pure.edges)[1]) {
    if (all(amat[pure.edges[i, 1], pure.edges[i, 2]] == amat[pure.edges[i,
                                                                        2], pure.edges[i, 1]])) {
      res$type[i] <- 0
      res$head[i] <- NA
      res$tail[i] <- NA
    } else {
      res$head[i] <- pure.edges[i, 2]
      res$tail[i] <- pure.edges[i, 1]
    }
  }
  s.pure.edges <- t(apply(pure.edges, 1, sort))
  ii <- order(s.pure.edges[, 1], s.pure.edges[, 2])
  res <- res[ii, ]
  res$xmin <- s.pure.edges[ii, 1]
  res$xmax <- s.pure.edges[ii, 2]
  return(res)
}
