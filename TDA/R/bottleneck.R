bottleneck <-
function(Diag1, Diag2, dimension = 1) {

  if (((class(Diag1) != "diagram" && class(Diag1) != "matrix" &&
      !is.data.frame(Diag1)) || NCOL(Diag1) != 3) &&
      (!is.numeric(Diag1) || length(Diag1) != 3)) {
    stop("Diag1 should be a diagram or a P by 3 matrix")
  }
  if (((class(Diag2) != "diagram" && class(Diag2) != "matrix" &&
      !is.data.frame(Diag2)) || NCOL(Diag2) != 3) &&
      (!is.numeric(Diag2) || length(Diag2) != 3)) {
    stop("Diag2 should be a diagram or a P by 3 matrix")
  }
  if (!is.numeric(dimension) || any(dimension < 0)) {
    stop("dimension should be a nonnegative integer or a vector of nonnegative integer")
  }

  if (is.numeric(Diag1)) {
    Diag1 <- matrix(Diag1, ncol = 3, dimnames = list(NULL, names(Diag1)))
  }
  if (is.numeric(Diag2)) {
    Diag2 <- matrix(Diag2, ncol = 3, dimnames = list(NULL, names(Diag2)))
  }

  bottleDistance <- rep(0, length(dimension))
  for (dimIdx in seq(along = dimension)) {
    Diag1Dim <- Diag1[which(Diag1[, 1] == dimension[dimIdx]), 2:3, drop = FALSE]
    Diag2Dim <- Diag2[which(Diag2[, 1] == dimension[dimIdx]), 2:3, drop = FALSE]

    if (length(Diag1Dim) == 0) {
      Diag1Dim <- array(c(dimension[dimIdx], 0, 0), dim = c(1, 3))
    }
    if (length(Diag2Dim) == 0) {
      Diag2Dim <- array(c(dimension[dimIdx], 0, 0), dim = c(1, 3))
    }

    bottleDistance[dimIdx] <- Bottleneck(Diag1Dim, Diag2Dim)
  }

  return (max(bottleDistance))
}
