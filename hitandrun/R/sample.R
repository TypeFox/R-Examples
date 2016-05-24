findFace <- function(x, constr) {
  stopifnot(length(x) == ncol(constr$constr))
  d <- constr$constr %*% x - constr$rhs
  which.max(d)
}

checkPolytope <- function(x0, constr, homogeneous, transform) {
  n <- if (!is.null(x0)) length(x0) else ncol(constr$constr)
  m <- nrow(constr$constr)

  # Verify preconditions
  stopifnot(n > homogeneous)
  stopifnot(n == ncol(constr$constr))
  stopifnot(m == length(constr$rhs))
  stopifnot(constr$dir == "<=")

  if (homogeneous) { # Change to homogeneous coordinates
    stopifnot(x0[n] == 1.0)
    list(n = n - 1,
         m = m,
         x0 = if (is.null(x0)) x0 else x0[1:(n - 1)],
         constr = list(constr = constr$constr[ , 1:(n - 1), drop=FALSE],
                       rhs = constr$rhs - constr$constr[ , n, drop=TRUE],
                       dir = constr$dir),
         transform = function(samples) {
           if (!is.null(transform)) {
             mat <- samples %*% t(transform[ , 1:(n - 1), drop=FALSE])
             t(t(mat) + transform[ , n, drop=TRUE])
           } else {
             cbind(samples, 1)
           }
         },
         xN = function(samples) {
           c(samples[nrow(samples), , drop=TRUE], 1)
         })
  } else {
    list(n = n,
         m = m,
         x0 = x0,
         constr = constr,
         transform = function(samples) {
           if (!is.null(transform)) {
             samples %*% t(transform)
           } else {
             samples
           }
         },
         xN = function(samples) {
           samples[nrow(samples), , drop=TRUE]
         })
  }
}

har <- function(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
  stopifnot(N %% thin == 0)
  args <- checkPolytope(x0, constr, homogeneous, transform)

  rval <- .Call("hitandrun_har", args$x0, args$constr$constr, args$constr$rhs, N, thin)

  list(samples=args$transform(rval),
       xN=args$xN(rval))
}

sab <- function(x0, i0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
  stopifnot(N %% thin == 0)
  args <- checkPolytope(x0, constr, homogeneous, transform)

  constr <- args$constr
  # normalize the constraints (required for shake-and-bake)
  for (i in 1:args$m) {
    norm <- sqrt(sum(constr$constr[i,]^2))
    constr$constr[i,] <- constr$constr[i,] / norm
    constr$rhs[i] <- constr$rhs[i] / norm
  }

  rval <- .Call("hitandrun_sab", args$x0, i0, constr$constr, constr$rhs, N, thin)

  list(samples=args$transform(rval[[1]]),
       xN=args$xN(rval[[1]]),
       faces=rval[[2]],
       iN=rval[[2]][length(rval[[2]])])
}

bbReject <- function(lb, ub, constr, N, homogeneous=FALSE, transform=NULL) {
  args <- checkPolytope(NULL, constr, homogeneous, transform)
  stopifnot(args$n == length(lb))
  stopifnot(args$n == length(ub))

  rval <- .Call("hitandrun_bbReject", lb, ub, args$constr$constr, args$constr$rhs, N)
  list(samples=args$transform(rval[[1]]),
       rejectionRate=rval[[2]])
}

simplex.sample <- function(n, N, sort=FALSE) {
  samples <- .Call("hitandrun_simplexSample", n, sort, N);
  list(samples=samples)
}

hypersphere.sample <- function(n, N) { 
  .Call("hitandrun_hypersphereSample", n, N)
}
