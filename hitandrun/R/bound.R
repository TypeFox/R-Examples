homogeneousCoordinateConstraint <- function(n) {
  list(constr=c(rep(0, n), 1), rhs=c(1), dir=c("="))
}

findInteriorPoint <- function(constr, homogeneous=FALSE, randomize=FALSE) {
  # max d, subject to
  # Ax + e = b ; slack on original constraints
  # d <= f_i e_i   ; minimum slack
  # where f_i = 1 if randomize = FALSE and  f_i ~ U(0.01, 1) otherwise
  n <- dim(constr$constr)[2] # number of variables
  m <- dim(constr$constr)[1] # number of constraints

  if (n == homogeneous) { # detect degenerate case
    return(rep(1, homogeneous))
  } else if (m < n) {
    stop("Underconstrained problem")
  }

  if (any(constr$dir != "<=")) {
    stop("findInteriorPoint only allows <= constraints")
  }

  f <- if (randomize) {
    runif(m, min=0.01, max=1)
  } else {
    rep(1, m)
  }
  ineq.constr <- matrix(0, nrow=m, ncol=n+m+1)
  ineq.constr[1:m,(n+1):(n+m)] <- -diag(f)
  ineq.constr[1:m,n+m+1] <- 1
  ineq.rhs <- rep(0, m)

  eq.constr <- matrix(0, nrow=m+homogeneous, ncol=n+m+1)
  eq.constr[1:m,1:n] <- constr$constr
  eq.constr[1:m,(n+1):(n+m)] <- diag(m)
  eq.rhs <- constr$rhs
  if (homogeneous == TRUE) {
    hom <- homogeneousCoordinateConstraint(n - 1)
    eq.constr[m+1,1:n] <- hom$constr
    eq.rhs <- c(eq.rhs, hom$rhs)
  }

  h <- rcdd::makeH(ineq.constr, ineq.rhs, eq.constr, eq.rhs)

  obj <- c(rep(0, n+m), 1) # maximize minimum slack
  sol <- rcdd::lpcdd(rcdd::d2q(h), rcdd::d2q(obj), minimize=FALSE)
  if (sol$solution.type != "Optimal") {
    stop(paste("No solution:", sol$solution.type))
  }
  sol <- rcdd::q2d(sol$primal.solution)
  if (sol[n+m+1] <= 0) stop("hitandrun::findInteriorPoint: infeasible constraints")
  sol[1:n]
}

findExtremePoints <- function(constr, homogeneous=FALSE) {
  n <- dim(constr$constr)[2]
  nh <- n

  if (any(constr$dir != "<=")) {
    stop("findInteriorPoint only allows <= constraints")
  }

  h <- if (homogeneous == TRUE) {
    n <- n - 1
    hom <- homogeneousCoordinateConstraint(n)
    rcdd::makeH(constr$constr, constr$rhs, hom$constr, hom$rhs)
  } else {
    rcdd::makeH(constr$constr, constr$rhs)
  }

  # for each of the n dimensions, solve 2 LPs to find the min/max
  findExtreme <- function(minimize) {
    function(i) {
      obj <- rep(0, nh)
      obj[i] <- 1
      sol <- rcdd::lpcdd(rcdd::d2q(h), rcdd::d2q(obj), minimize=minimize)
      if (sol$solution.type != "Optimal") {
        stop(paste("No solution:", sol$solution.type))
      }
      rcdd::q2d(sol$primal.solution)
    }
  }
  t(cbind(sapply(1:n, findExtreme(TRUE)), sapply(1:n, findExtreme(FALSE))))
}

findVertices <- function(constr, homogeneous=FALSE) {
  if (any(constr$dir != "<=")) {
    stop("findInteriorPoint only allows <= constraints")
  }

  h <- if (homogeneous == TRUE) {
    n <- dim(constr$constr)[2]
    hom <- homogeneousCoordinateConstraint(n - 1)
    rcdd::makeH(constr$constr, constr$rhs, hom$constr, hom$rhs)
  } else {
    rcdd::makeH(constr$constr, constr$rhs)
  }

  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  # Check that the output are vertices, not other things that would indicate a problem
  if (v[,1] != "0" || v[,2] != "1") {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }

  # Return the vertices only 
  v[ , -c(1,2), drop=FALSE]
}

# Generate seed point from constraints
createSeedPoint <- function(constr, homogeneous=FALSE, randomize=FALSE,
    method="slacklp") {
  stopifnot(method %in% c("slacklp", "vertices"))

  p <- if (method == "slacklp") {
    findInteriorPoint(constr, homogeneous, randomize)
  } else { # method == "vertices"
    n <- dim(constr$constr)[2]
    if (homogeneous == TRUE) {
      n <- n - 1
    }

    extreme <- findVertices(constr, homogeneous)

    # starting point 
    m <- nrow(dim(extreme))
    if (randomize == TRUE) { # random weighting
      w <- as.vector(simplex.sample(m, 1)$samples)
      apply(extreme, 2, function(row) { sum(w * row) })
    } else { # mean: approximation of centroid
      apply(extreme, 2, mean)
    }
  }

  if (homogeneous == TRUE) {
    p[length(p)] <- 1.0 # eliminate floating point imprecision
  }
  p
}

# Create a bounding box given constraints in (n-1)
createBoundBox <- function(constr, homogeneous=FALSE) {
  n <- dim(constr$constr)[2]
  extreme <- findExtremePoints(constr, homogeneous)
  if (homogeneous == TRUE) {
    extreme <- extreme[ , 1:(n-1), drop=FALSE]
  }
  # upper and lower bounds for each dimension in the (n-1) basis
  lb <- apply(extreme, 2, min)
  ub <- apply(extreme, 2, max)
  list(lb=lb, ub=ub)
}
