# Morris' simplex sub-routines. (See main file morris.R)
#
# Gilles Pujol 2007-2008
# Modified by Frank Weber (2016): support model functions
# returning a matrix or a 3-dimensional array.

simplex.reg <- function(p) {
  # generates the matrix of a regular simplex, of edge length = 1,
  # centered on the origin
  S <- matrix(0, nrow = p + 1, ncol = p)
  S[2,1] <- 1
  for (i in 3 : (p + 1)) {
    for (j in 1 : (i - 2)) {
      S[i,j] <- mean(S[1 : (i - 1), j])
    }
    S[i, i - 1] <- sqrt(1 - sum(S[i, 1 : (i - 2)]^2))
  }
  scale(S, scale = FALSE)
}


simplex.rect <- function(p) {
  # generates the matrix of a rectangle ("orthonormal") simplex
  S <- matrix(0, nrow = p + 1, ncol = p)
  for (i in 1:p) {
    S[i+1,i] <- 1
  }
  S
}


plane.rot <- function(p, i, j, theta) {
  # matrix of the plane (i,j)-rotation of angle theta in dimension p
  R <- diag(nrow = p)
  R[c(i,j), c(i,j)] <- matrix(c(cos(theta), sin(theta), - sin(theta), cos(theta)), nrow = 2)
  return(R) 
}

random.simplexes <- function(p, r, min = rep(0, p), max = rep(1, p), h = 0.25) {
  # generates r random simplexes
  
  # Check if p >= 2:
  if(p < 2){
    stop("At least 2 factors have to be analyzed.")
  }
  X <- matrix(nrow = r * (p + 1), ncol = p)
  
  # initial random rotation matrix
  R <- diag(nrow = p, ncol = p)
  ind <- combn(p, 2)
  theta <- runif(choose(p, 2), min = 0, max = 2 * pi)
  for (i in 1 : choose(p, 2)) {
    R <- plane.rot(p, ind[1,i], ind[2,i], theta[i]) %*% R
  }
  
  # reference simplex
  S.ref <- R %*% (h * t(simplex.reg(p)))
  
  for (i in 1 : r) {
    # one more random plane rotation of the reference simplex
    ind <- sample(p, 2)
    theta <- runif(1, min = 0, max = 2 * pi)
    S.ref <- plane.rot(p, ind[1], ind[2], theta) %*% S.ref
    
    # translation to put the simplex into the domain
    tau.min <- min - apply(S.ref, 1, min)
    tau.max <- max - apply(S.ref, 1, max)
    X[ind.rep(i, p),] <- t(S.ref + runif(n = p, min = tau.min, max = tau.max))
  }
  
  return(X)
}


ee.simplex <- function(X, y) {
  # compute the elementary effects for a simplex design
  p <- ncol(X)
  r <- nrow(X) / (p + 1)
  if(class(y) == "numeric"){
    one_i_vector <- function(i){
      j <- ind.rep(i, p)
      return(solve(cbind(as.matrix(X[j, ]), rep(1, p + 1)), y[j])[1:p])
    }
    ee <- vapply(1:r, one_i_vector, FUN.VALUE = numeric(p))
    ee <- t(ee)
    # "ee" is now a (r times p)-matrix.
  } else if(class(y) == "matrix"){
    one_i_matrix <- function(i){
      j <- ind.rep(i, p)
      return(solve(cbind(as.matrix(X[j, ]), rep(1, p + 1)), 
                   y[j, , drop = FALSE])[1:p, , drop = FALSE])
    }
    ee <- vapply(1:r, one_i_matrix, 
                 FUN.VALUE = matrix(0, nrow = p, ncol = ncol(y)))
    # Transpose "ee" (an array of dimensions c(p, ncol(y), r)) to an array of
    # dimensions c(r, p, ncol(y)) (for better consistency with the standard 
    # case that "class(y) == "numeric""):
    ee <- aperm(ee, perm = c(3, 1, 2))
  } else if(class(y) == "array"){
    one_i_array <- function(i){
      j <- ind.rep(i, p)
      ee_per_3rd_dim <- sapply(1:(dim(y)[3]), function(idx_3rd_dim){
        y_j_matrix <- y[j, , idx_3rd_dim]
        # Correction needed if "dim(y)[2] == 1", so "y_j_matrix" has been 
        # dropped to a vector:
        if(class(y_j_matrix) == "numeric"){
          y_j_matrix <- matrix(y_j_matrix)
        }
        # Here, the result of "solve(...)" is a (p times dim(y)[2])-matrix:
        solve(cbind(as.matrix(X[j, ]), rep(1, p + 1)), 
              y_j_matrix)[1:p, , drop = FALSE]
      }, simplify = "array")
      # "ee_per_3rd_dim" is an array of dimensions c(p, dim(y)[2], dim(y)[3]).
      # Assign the corresponding names for the third dimension:
      dimnames(ee_per_3rd_dim)[[3]] <- dimnames(y)[[3]]
      return(ee_per_3rd_dim)
    }
    ee <- sapply(1:r, one_i_array, simplify = "array")
    # "ee" is an array of dimensions c(p, dim(y)[2], dim(y)[3], r), so it is
    # transposed to an array of dimensions c(r, p, dim(y)[2], dim(y)[3]):
    ee <- aperm(ee, perm = c(4, 1, 2, 3))
  }
  return(ee)
}
