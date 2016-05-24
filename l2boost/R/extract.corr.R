# This is a hidden function of the l2boost package.

# extract correlations: uses a smart reduction for elastic net modification

# @param x design matrix (possibly augmented for elasticBoosting)
# @param l currently selected direction
# @param enet elasticBoost indicator, required since we only scale the non-augmented design matrix
# @param n.org length of the non-augmented design matrix

extract.corr <- function(x, l, enet, n.org) {
  if (enet) {
    e.c <- c(t(x[1:n.org, ]) %*% x[1:n.org, l])
    e.c[l] <- 1
    e.c
  }
  else {
    c(t(x) %*% x[, l])
  }
}
