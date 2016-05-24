# Space-filling optimization of a Morris design, either OAT or
# simplex-based. (See main file morris.R)
#
# Gilles Pujol 2007


haussdorf.distance <- function(x, set1, set2) {
# x: matrix of points.
# set1: indices of points (in x) of the first group.
# set2: indices of points (in x) of the second group.
# returns: the Haussdorf distance between the two sets of points.
  n1 <- length(set1)
  n2 <- length(set2)
  d <- matrix(nrow = n1, ncol = n2)
  for (i1 in 1 : n1) {
    for (i2 in 1 : n2) {
      d[i1,i2] <- sqrt(sum((x[set1[i1],] - x[set2[i2],])^2))
    }
  }
  return(max(mean(apply(d, 1, min)), mean(apply(d, 2, min))))
}


kennard.stone <- function(dist.matrix, n) {
# Kennard & Stone algorithm (1969).
# dist.matrix: distance matrix (N * N) (cf help(dist)).
# n: number of points to keep (n < N).
# returns: the indices of the n chosen points.
  out <- numeric(n)
  out[1] <- 1
  for (i in 2 : n) {
    tmp <- dist.matrix[out, -out, drop = FALSE]
    # Remark: drop = FALSE since 'out' is of length 1 at the first
    # iteration, cf help(Extract) for the meaning of 'drop'
    out[i] <- (1 : nrow(dist.matrix))[-out][which.max(apply(tmp, 2, min))]
  }
  return(out)
}


morris.maximin <- function(x, r) {
# Select r repetitions (out of the R ones of the "morris" design x)
# that are "space-filling".
# returns: the indices (in 1:R) of the r selected repetitions.
  p <- ncol(x)
  R <- nrow(x) / (p + 1)
  d <- matrix(0, nrow = R, ncol = R)
  for (i in 1 : (R - 1)) {
    for (j in (i + 1) : R) {
      d[i,j] <- d[j,i] <- haussdorf.distance(x, ind.rep(i, p), ind.rep(j, p))
    }
  }
  kennard.stone(d, r)
}
