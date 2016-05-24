
# Produces a matrix of distances between points.

# Not exported, no sanity checks, no documentation.

# Arguments:
# A, B : 2-column MATRICES giving the coordinates of each set of points.

# Returns:
# A matrix of distances # A in rows, B in columns? Named?

# This is Andy Royle's original function:
# (Mike changed x/y to A/B to avoid confusion with x=east, y=north)

# `e2dist` <- function (A, B) {
  # i <- sort(rep(1:nrow(B), nrow(A)))
  # dvec <- sqrt((A[, 1] - B[i, 1])^2 + (A[, 2] - B[i, 2])^2)
  # matrix(dvec, nrow = nrow(A), ncol = nrow(B), byrow = FALSE)
# }

# This is Mike Meredith's version, which preserves row names
#   and may run faster:

e2dist <- function(A, B)  {
  xdif <- outer(A[, 1], B[, 1], "-")
  ydif <- outer(A[, 2], B[, 2], "-")
  sqrt(xdif^2 + ydif^2)
}

