offdiagonal <- function(g, deg = NULL) {

  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(deg))
    deg <- graph::degree(g)
  n <- numNodes(g)
  M <- adjacencyMatrix(g)

  nr <- max(deg) + 1
  c_nodecor <- matrix(0, nrow=nr, ncol=nr)
  for (l in 1:n)
    for (p in 1:n)
      if (M[l, p] != 0 && deg[p] >= deg[l])
        c_nodecor[deg[l]+1, deg[p]+1] <- c_nodecor[deg[l]+1, deg[p]+1] + 1

  # sum over upper-triangle "parallels" to the main diagonal
  a_offdiag <- sapply(0:(nr - 1), function(delta) {
    sum(c_nodecor[cbind(1:(nr-delta), (delta+1):nr)])
  })

  prob_a <- a_offdiag / sum(a_offdiag)
  sum_prob <- sum(sapply(prob_a, function(p) if (p == 0) 0 else p * log(p)))

  -sum_prob / log(n-1)
}
