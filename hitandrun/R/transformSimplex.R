# calculate the Moore-Penrose pseudo-inverse of matrix
pseudoinverse <- function (m) {
  msvd <- svd(m)
  if (length(msvd$d) == 0) {
    array(0, dim(m)[2:1])
  } else {
    s <- matrix(0, nrow=length(msvd$d), ncol=length(msvd$d))
    diag(s)[msvd$d != 0] <- 1/msvd$d[msvd$d != 0]
    msvd$v %*% (s %*% t(msvd$u))
  }
}

# generate basis for solution space of equality constraints
solution.basis <- function(constr) {
  stopifnot(all(constr$dir == '='))
  A <- constr$constr
  b <- constr$rhs
  n <- ncol(A)

  A.inv <- pseudoinverse(A)

  if (isTRUE(all.equal(diag(n), A.inv %*% A))) {
    list(basis=matrix(0, nrow=ncol(A), ncol=0),
         translate=as.vector(A.inv %*% b))
  } else {
    the.qr <- qr(diag(n) - A.inv %*% A)
    list(basis=qr.Q(the.qr)[, 1:the.qr$rank, drop=FALSE],
         translate=as.vector(A.inv %*% b))
  }
}

# create basis for (translated) n-dim simplex
simplex.basis <- function(n) {
  b <- rbind(diag(n-1), rep(-1, n-1))
  list(basis=qr.Q(qr(b)),
       translate=rep(1/n, n))
}

# Generate a projection matrix that transforms an (n-1) dimensional vector in
# homogeneous coordinate representation to an n-dimensional weight vector.
createTransform <- function(basis, inverse=FALSE, keepHomogeneous=inverse) {
  # add one extra element to vectors in each basis (homogeneous coordinate
  # representation)
  translate <- if (inverse == FALSE) basis$translate else -basis$translate
  n <- length(translate)
  basis <- basis$basis
  basis <- rbind(cbind(basis, rep(0, n)), c(rep(0, ncol(basis)), 1))

  # create translation matrix (using homogenous coordinates)
  translation <- rbind(cbind(diag(n), translate), c(rep(0, n), 1))

  # homogeneous coordinate elimination
  nh <- if (inverse == FALSE) { nrow(basis) } else { ncol(basis) }
  elimHom <- if (keepHomogeneous) {
    diag(nh)
  } else {
    cbind(diag(nh - 1), rep(0, nh - 1))
  }

  transform <- if (inverse == FALSE) {
    # successively apply basis transformation and translation
    elimHom %*% translation %*% basis
  } else {
    # successively apply translation and basis transformation
    elimHom %*% t(basis) %*% translation
  }

  unname(transform)
}

# Generate a projection matrix that transforms an (n-1) dimensional vector in
# homogeneous coordinate representation to an n-dimensional weight vector.
simplex.createTransform <- function(n, inverse=FALSE, keepHomogeneous=inverse) {
  createTransform(simplex.basis(n), inverse, keepHomogeneous)
}

transformConstraints <- function(transform, constr) {
  list(
    constr = constr$constr %*% transform,
    dir = constr$dir,
    rhs = constr$rhs)
}

# translate the n-dimensional constraints to the (n-1)-dimensional space
# transform: transform created by simplex.createTransform
# userConstr: additional constraints
simplex.createConstraints <- function(transform, userConstr=NULL) {
  n <- nrow(transform)

  # basic constraints defining the (n-1)-dimensional simplex
  constr <- list(
    constr = diag(rep(-1, n)), # -1*w[i] <= 0
    dir = rep('<=', n),
    rhs = rep(0, n))

  # user constraints
  if (!is.null(userConstr)) {
    stopifnot(userConstr$dir == "<=")
    constr <- mergeConstraints(constr, userConstr)
  }

  transformConstraints(transform, constr)
}
