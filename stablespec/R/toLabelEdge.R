toLabelEdge <- function(adjMat, consMatrix) {
  # all credit goes to Markus Kalisch (kalisch@stat.math.ethz.ch) and Diego Colombo who
  # originally implemented this. We borrow this function and
  # slightly modified so as to meet our need.
  #
  # Kalisch, M., Machler, M., Colombo, D., Maathuis, M. H., &
  # Buehlmann, P. (2012). Causal inference using graphical models
  # with the R package pcalg. Journal of Statistical Software, 47(11), 1-26.

  edge.df <- orderEdge(adjMat)
  lab <- rep(NA, dim(edge.df)[1])
  edge.df <- edge.df[order(edge.df$order), ]
  Head <- edge.df$head
  Tail <- edge.df$tail

  for (i in 1:nrow(consMatrix)) {
    consTail <- consMatrix[i, 1]
    consHead <- consMatrix[i, 2]
    ind <- which(edge.df$tail == consTail & edge.df$head == consHead)
    lab[ind] <- TRUE
    lab[which(lab != TRUE)] <- NA
  }

  while (any(ina <- is.na(lab))) {
    x.y <- which(ina)[1]
    x <- Tail[x.y]
    y <- Head[x.y]
    y.is.head <- Head == y
    e1 <- which(Head == x & lab)
    for (ee in e1) {
      w <- Tail[ee]
      if (any(wt.yh <- w == Tail & y.is.head))
        lab[wt.yh] <- TRUE
      else {
        lab[y.is.head] <- TRUE
        break
      }
    }

    cand <- which(y.is.head & Tail != x)
    if (length(cand) > 0) {
      valid.cand <- rep(FALSE, length(cand))
      for (iz in seq_along(cand)) {
        z <- Tail[cand[iz]]
        if (!any(Tail == z & Head == x))
          valid.cand[iz] <- TRUE
      }
      cand <- cand[valid.cand]
    }
    lab[which(y.is.head & is.na(lab))] <- (length(cand) > 0) # label with FALSE or "REVERSIBLE"
  }

  edge.df$label <- lab
  return(edge.df)
}

