library(pcalg)

# Perform tests with only two  DAGs, but with permutations of the vertices
# (to check for a bug present till pcalg-2.0.8)

n.perm <- 5

set.seed(123)

# A: adjacency matrix of DAG;
# B: adjacency matrix of CPDAG
# Setting 3 by courtesy of Jonas Peters: in pcalg <= 2.0.8,
# setting i = 3, k = 3 failed
A <- list(
  matrix(c(0,0,0,0,1, 0,0,1,0,1, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,1,0,0,0, 0,0,0,1,0, 0,0,0,1,0, 0,0,0,0,1, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,0,0,0, 1,0,0,0, 1,1,0,0, 1,1,1,0), 4, 4))
B <- list(
  matrix(c(0,0,0,0,1, 0,0,1,0,1, 0,1,0,1,0, 0,0,1,0,0, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,1,0,0,0, 1,0,0,1,0, 0,0,0,1,0, 0,0,0,0,1, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,1,1,1, 1,0,1,1, 1,1,0,1, 1,1,1,0), 4, 4))

for (i in 1:length(A)) {
  for (k in 1:n.perm) {
    p <- nrow(A[[i]])
    
    ind <- if(k == 1) 1:p else sample.int(p)
    
    g <- as(A[[i]][ind, ind], "graphNEL")
    pdag <- dag2cpdag(g)
    B.hat <- as(pdag, "matrix")
    if (!all(B.hat == B[[i]][ind, ind])) {
      stop(sprintf("True CPDAG not found! (setting: i = %d, k = %d)", i, k))
    }
    
    # par(mfrow = c(1, 2))
    # plot(g)
    # plot(pdag)
  }
}
