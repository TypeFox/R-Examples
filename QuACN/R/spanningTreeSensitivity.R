spanningTreeSensitivity <- function(g, one.eds=NULL) {

  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(one.eds))
    one.eds <- edgeDeletedSubgraphs(g)

  n <- numNodes(g)
  m <- numEdges(g)

  # number of spanning trees in g
  lap <- laplaceMatrix(g)
  nST_g <- det(lap[2:n, 2:n])

  # number of spanning trees in each subgraph
  nST_1e <- sapply(one.eds, function(M_1e) {
    diag_1e <- diag(rowSums(M_1e, na.rm = FALSE, dims = 1))
    lap_1e <- diag_1e - M_1e
    det(lap_1e[2:n, 2:n])
  })

  m_cu <- n^1.68 - 10

  # spanning tree sensitivity complexity measure
  sens <- nST_g - nST_1e
  sens <- sort(unique(sens))
  sens <- sens[sens != 0]
  sens <- sens - min(sens) + 1
  prob <- sens / sum(sens)
  STS <- -sum(prob * log(prob)) / log(m_cu)

  # STSD complexity measure
  sens_diff <- unique(diff(sens))
  prob <- sens_diff / sum(sens_diff)
  prob <- prob[prob != 0]
  STSD <- -sum(prob * log(prob)) / log(m_cu)

  list(`STS` = STS, `STSD` = STSD)
}
