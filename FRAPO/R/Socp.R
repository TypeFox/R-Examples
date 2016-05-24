##
## Functions for solving SOCP
##
Socp <- function (f, A, b, C, d, N, x = NULL, z = NULL, w = NULL, control = list()){
  cpos <- SocpControl()
  if (length(control)){
    nms <- names(control)
    if (!is.list(control) || is.null(nms))
      stop("control argument must be a named list")
    pos <- pmatch(nms, names(cpos))
    if (any(nap <- is.na(pos))){
      warning(paste("unrecognized control element(s) named `", paste(nms[nap], collapse = ", "), "' ignored", sep = ""))
      pos <- pos[!nap]
      control <- control[!nap]
    }
    if (length(control))
      control <- c(control, cpos[-pos])
    else control <- cpos
    } else {
      control <- cpos
    }
  abs.tol <- as.double(control$abs.tol)
  rel.tol <- as.double(control$rel.tol)
  target <- as.double(control$target)
  max.iter <- as.integer(control$max.iter)
  Nu <- as.double(control$Nu)
  out.mode <- as.integer(control$out.mode)
  BigM.K <- as.integer(control$BigM.K)
  BigM.iter <- as.integer(control$BigM.iter)
  stopifnot(Nu >= 1)
  stopifnot(out.mode %in% 0:3)
  L <- length(N)
  n <- length(f)
  iter <- max.iter
  mhist <- 0
  nhist <- 0
  ndbl <- 0
  nint <- 0
  info <- 0
  if (!is.null(d)){
    stopifnot(any(dim(A) == c(sum(N), n)))
    stopifnot(length(b) == sum(N))
    stopifnot(any(dim(C) == c(L, n)))
    stopifnot(length(d) == L)
    if (!is.null(x)) stopifnot(length(x) == n)
    if (!is.null(z)) stopifnot(length(z) == sum(N))
    if (!is.null(z)) stopifnot(length(w) == L)
    Ni <- seq(from = 1, length = N[1])
    AA <- rbind(A[Ni, ], C[1, ])
    bb <- c(b[Ni], d[1])
    if (!is.null(z))
      zz <- c(z[Ni], w[1])
    if(L > 1){
      Ncum <- cumsum(N) + 1
      for (i in seq(L, from = 2)){
        Ni <- seq(from = Ncum[i-1], length = N[i])
        AA <- rbind(AA, A[Ni, ], C[i, ])
        bb <- c(bb, b[Ni], d[i])
        if (!is.null(z)) zz <- c(zz, z[Ni], w[i])
      }
    }
    A <- AA
    b <- bb
    if (!is.null(z)) z <- zz
    N <- N + 1
  }
  if (is.null(x)) x <- SocpPhase1(f, A, b, N, control)
  if (is.null(z)) z <- SocpPhase2(f, A, b, N, x, control)
  socp.getwork <- .C("socp_getwork", as.integer(L), as.integer(N), as.integer(n), as.integer(max.iter), as.integer(out.mode), as.integer(mhist), as.integer(nhist), as.integer(ndbl), as.integer(nint), PACKAGE = "FRAPO")
  mhist <- socp.getwork[[6]]
  nhist <- socp.getwork[[7]]
  ndbl <- socp.getwork[[8]]
  nint <- socp.getwork[[9]]
  hist <- rep(0, nhist * mhist + 1)
  dblwork <- rep(0, ndbl + 1)
  intwork <- rep(0, nint + 1)
  socp <- .C("socp", as.integer(L), as.integer(N), as.integer(n), as.double(f), as.double(as.vector(A)), as.double(b), as.double(x), as.double(z), as.double(abs.tol), as.double(rel.tol), as.double(target), as.integer(iter), as.double(Nu), as.integer(info), as.integer(out.mode), as.double(hist), as.double(dblwork), as.integer(intwork), PACKAGE = "FRAPO")
  ans <- list(x = socp[[7]])
  if (!is.null(d)){
    idx <- cumsum(N)
    ans$z <- socp[[8]][-idx]
    ans$w <- socp[[8]][idx]
  } else {
    ans$z <- socp[[8]]
  }
  ans$iter <- socp[[12]]
  ans$hist <- socp[[16]]
  ans$convergence <- as.logical(if (socp[[14]] %in% 1:3) TRUE else FALSE)
  ans$info <- socp[[14]]
  ans$message <- switch(as.character(socp[[14]]), `0` = "Error (0)", `1` = "Absolute tolerance (1)", `2` = "Relative convergence (2)", `3` = "Target value (achieved or unachievable) (3)", `4` = "Maximum iterations (4)")
  if (ans$info == 1 | ans$info == 2){
    e2 <- sum((t(A) %*% socp[[8]] - f)^2)^0.5/(sum(socp[[8]]^2))^0.5
    if (e2 > 1e-06) warning("Dual is not in the feasible hyperplane, solution may be wrong.")
  }
  return(ans)
}
##
## Initialising parameters
##
SocpPhase1 <- function (f, A, b, N, control){
  abs.tol <- as.double(control$abs.tol)
  rel.tol <- as.double(control$rel.tol)
  target <- as.double(control$target)
  max.iter <- as.integer(control$max.iter)
  Nu <- as.double(control$Nu)
  out.mode <- as.integer(control$out.mode)
  BigM.K <- as.integer(control$BigM.K)
  BigM.iter <- as.integer(control$BigM.iter)
  n <- length(f)
  m <- nrow(A)
  L <- length(N)
  SU <- matrix(0, nrow = m, ncol = L)
  ST <- matrix(0, nrow = m, ncol = L)
  T <- rep(0, m)
  k <- 1
  for (i in seq(L)){
    if (N[i] > 1)
      SU[k:(k + N[i] - 2), i] = matrix(1, nrow = N[i] - 1, ncol = 1)
    ST[k + N[i] - 1, i] <- 1
    T[k + N[i] - 1] <- 1
    k <- k + N[i]
  }
  slacka <- max(rowSums(abs(cbind(A, b))))
  x <- rep(0, n)
  u0 <- A %*% x + b
  k <- -t(ST) %*% u0 + sqrt(t(SU) %*% (u0^2))
  alpha <- slacka + max(k)
  x0 <- qr.solve(A, T)
  if (sum(abs(T - A %*% x0)^2)^(1/2) < 10 * .Machine$double.eps * sqrt(m)){
    x <- alpha * x0
  } else {
    if (alpha >= 0) {
      x1 <- c(x, alpha)
      f1 <- c(rep(0, n), 1)
      A1 <- cbind(A, T)
      fit <- Socp(f1, A1, b, NULL, NULL, N, x1, NULL, NULL, control = list(target = 0, rel.tol = -1))
      if (fit$x[n + 1] >= 0)
        stop("Phase 1 failed, alpha>=0")
      x <- fit$x[1:n]
    }
  }
  return(x)
}
##
## Initialising values of z in dual form
##
SocpPhase2 <- function (f, A, b, N, x, control){
  abs.tol <- as.double(control$abs.tol)
  rel.tol <- as.double(control$rel.tol)
  target <- as.double(control$target)
  max.iter <- as.integer(control$max.iter)
  Nu <- as.double(control$Nu)
  out.mode <- as.integer(control$out.mode)
  BigM.K <- as.integer(control$BigM.K)
  BigM.iter <- as.integer(control$BigM.iter)
  m <- nrow(A)
  L <- length(N)
  SU <- matrix(0, nrow = m, ncol = L)
  ST <- matrix(0, nrow = m, ncol = L)
  T <- rep(0, m)
  k <- 1
  for (i in seq(L)){
    if (N[i] > 1) SU[k:(k + N[i] - 2), i] = rep(1, N[i] - 1)
    ST[k + N[i] - 1, i] <- 1
    T[k + N[i] - 1] <- 1
    k <- k + N[i]
  }
  BigM.repeat <- max.iter/BigM.iter
  normA <- max(colSums(abs(A)))
  normf <- max(abs(f))
  slackb <- 0.1 * normf/normA
  u <- A %*% x + b
  B0 <- BigM.K * (T %*% u)
  z0 <- qr.solve(t(A), f)
  k <- -t(ST) %*% z0 + sqrt(t(SU) %*% (z0^2))
  beta <- slackb + max(0, max(k))
  z = z0 + T * beta
  z1 = c(z, beta)
  A1 <- rbind(A, -T %*% A)
  b1 <- c(b, B0 - T %*% b)
  N1 <- c(N, 1)
  i <- 0
  info <- 4
  while (info == 4 & i < BigM.repeat){
    i <- i + 1
    fit <- Socp(f, A1, b1, NULL, NULL, N1, x, z1, NULL, control)
    info <- fit$info
    u <- A %*% fit$x + b
    t <- T %*% u
    if (BigM.K * t > B0){
      b1[m + 1] <- b1[m + 1] + BigM.K * t - B0
      B0 <- BigM.K * t
    }
  }
  e1 <- fit$z[m + 1]/((sum(fit$z[1:m]^2))^(1/2))
  if ((rel.tol >= 0) && (e1 > max(rel.tol, 1e-04)))
    warning("Warning: big M bound appears to be active, solution may be wrong.")
  z <- fit$z[1:m]
  return(z)
}
##
## controlling behaviour of algorithm
##
SocpControl <- function(abs.tol = 1e-18, rel.tol = 1e-16, target = 0, max.iter = 500, Nu = 10, out.mode = 0, BigM.K = 2, BigM.iter = 5){
  list(abs.tol = abs.tol, 
        rel.tol = rel.tol,
        target = target,
        max.iter = max.iter,
        Nu = Nu,
        out.mode = out.mode,
        BigM.K = BigM.K,
        BigM.iter = BigM.iter)
}
