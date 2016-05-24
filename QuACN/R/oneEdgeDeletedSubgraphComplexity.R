oneEdgeDeletedSubgraphComplexity <- function(g, one.eds=NULL) {


  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if(numEdges(g)==0)
    stop("No edges in current graph object")
  #################################
  
  if (is.null(one.eds))
    one.eds <- edgeDeletedSubgraphs(g)

  n <- numNodes(g)
  count <- length(one.eds)

  # number of spanning trees in g
  lap <- laplaceMatrix(g)
  nST_g <- det(lap[2:n, 2:n])

  # number of spanning trees in each subgraph and
  # eigenvalues of Laplacian and signless Laplacian of each subgraph
  data <- lapply(one.eds, function(M_1e) {
    diag_1e <- diag(rowSums(M_1e, na.rm = FALSE, dims = 1))
    lap_1e <- diag_1e - M_1e
    nST_1e <- det(lap_1e[2:n, 2:n])
    EV_lap_1e <- as.double(eigen(lap_1e, only.values = TRUE)$values)
    signless_lap_1e <- diag_1e + M_1e
    EV_signless_lap_1e <- as.double(eigen(signless_lap_1e, only.values = TRUE)$values)
    list(nST = nST_1e, EV_lap = EV_lap_1e, EV_signless_lap = EV_signless_lap_1e)
  })

  # number of subgraphs that are nonisomorphic in terms of
  # different spanning trees and spectrum of Laplacian and signless Laplacian
  sST <- 0
  sSpec <- 0
  for (k in 1:(count-1)) {
    for (l in (k+1):count) {
      if (data[[k]]$nST == data[[l]]$nST) {
        sST <- sST + 1
        break
      }
    }
    for (l in (k+1):count) {
      if (setequal(data[[k]]$EV_lap, data[[l]]$EV_lap) &&
          setequal(data[[k]]$EV_signless_lap, data[[l]]$EV_signless_lap)) {
        sSpec <- sSpec + 1
        break
      }
    }
  }
  N_1eST <- count - sST
  N_1eSpec <- count - sSpec

  # one-edge-deleted subgraph and spectrum complexity
  m_cu <- n^1.68 - 10
  C_1eST <- (N_1eST - 1) / (m_cu - 1)
  C_1eSpec <- (N_1eSpec - 1) / (m_cu - 1)

  list(`C_1eST` = C_1eST, `C_1eSpec` = C_1eSpec)
}
