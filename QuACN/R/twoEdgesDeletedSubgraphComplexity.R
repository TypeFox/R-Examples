twoEdgesDeletedSubgraphComplexity <- function(g, two.eds=NULL) {

  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(two.eds))
    two.eds <- edgeDeletedSubgraphs(edgeDeletedSubgraphs(g))

  n <- numNodes(g)
  count <- length(two.eds)

  data <- lapply(two.eds, function(M_2e) {
    diag_2e <- diag(rowSums(M_2e, na.rm = FALSE, dims = 1))
    lap_2e <- diag_2e - M_2e
    nST_2e <- det(lap_2e[2:n, 2:n])
    EV_lap_2e <- as.double(eigen(lap_2e, only.values = TRUE)$values)
    signless_lap_2e <- diag_2e + M_2e
    EV_signless_lap_2e <- as.double(eigen(signless_lap_2e, only.values = TRUE)$values)
    list(nST = nST_2e, EV_lap = EV_lap_2e, EV_signless_lap = EV_signless_lap_2e)
  })

  sSpec <- 0
  for (k in 1:(count-1)) {
    for (l in (k+1):count) {
      if (setequal(data[[k]]$EV_lap, data[[l]]$EV_lap) &&
          setequal(data[[k]]$EV_signless_lap, data[[l]]$EV_signless_lap)) {
        sSpec <- sSpec + 1
        break
      }
    }
  }

  N_2eSpec <- count - sSpec
  m_cu <- n^1.68 - 10
  C_2eSpec <- (N_2eSpec - 1) / (choose(m_cu, 2) - 1)

  C_2eSpec
}
