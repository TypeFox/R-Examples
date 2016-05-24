list.ancestors <- function (parent, child, node) {
    pvector <- numeric(max(parent))
    pvector[child] <- parent
    anc <- function(pvector, node) {
        res <- numeric(0)
        repeat {
            anc <- pvector[node]
            if (anc == 0) 
              break
            res <- c(res, anc)
            node <- anc
        }
        res
    }
    return(anc(pvector, node))
}

depths <- function (parent, child) {
  root   <- min(parent)
  depth  <- integer(max(parent))
  for (i in 1:length(parent)) {
    depth[child[i]] <- depth[parent[i]] + 1L
  }
  as.integer(depth)
}

tci <- function (tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nTip   <- length(tree$tip.label)
  depth  <- depths(parent, child)
  ancestors <- lapply(1:nTip, function(node) list.ancestors (parent, child, node))
  lca.depth <- vapply(1:nTip, function(i) {
    vapply(1:nTip, function (j) {
      anc.i <- ancestors[[i]]
      anc.j <- ancestors[[j]]
      lca <- max(anc.i[anc.i %in% anc.j])
      depth[lca]
    }, integer(1))
  }, integer(nTip))
  return (sum(lca.depth[upper.tri(lca.depth)]))
}

tci.context <- function (tree) {
  n <- length(tree$tip.label)
  maximum <- choose(n, 3)
  minimum <- mci(n)
  # Theorem 17
  uniform.expected <- (1 / 2) * choose(n, 2) * ((dfact((2 * n) - 2) / dfact((2 * n) - 3)) - 2)
  yule.expected    <- n * (n + 1) - (2 * n * H(n))
  yule.variance    <- ((1 / 12) * (n^4 - (10 * n^3) + (131 * n^2) - (2 * n))) - (4 * n^2 * H2(n)) - (6 * n * H(n))
  return (data.frame(maximum,   minimum,   uniform.expected,  yule.expected,  yule.variance))
}

H <- function (n) sum(1 / (1:n))
H2 <- function (n) sum(1 / (1:n)^2)

dfact <- function (n) {
  if (n < 2) return (1)
  n * dfact(n - 2)
}

mci <- function (n) { # Lemma 14 in Mir er al 2013
  if (n < 3) return (0)
  ceiling(mci(n/2)) + floor (mci(n/2)) + choose(ceiling(n/2), 2) + choose(floor(n/2), 2)
}

