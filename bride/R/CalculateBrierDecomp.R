#the main function
CalculateBrierDecomp <- function(p, y, n.bins) {

  n <- length(p)

  #p.binning (vector of length n with entries equal to bin-no. of p_n)
  p.breaks <- seq(0, 1, length.out = n.bins + 1)
  p.binning <- cut(p, breaks=p.breaks, include.lowest=TRUE, ordered_result=TRUE)
  p.binning <- as.numeric(p.binning)

  ##construct matrices and column sums
  m.ind <- matrix(c(seq(n), p.binning), nrow=n)
  m.A <- matrix(0, nrow=n, ncol=n.bins)
  m.A[m.ind] <- 1
  cs.A <- colSums(m.A)
  m.B <- matrix(0, nrow=n, ncol=n.bins)
  m.B[m.ind] <- y
  cs.B <- colSums(m.B)
  m.C <- matrix(0, nrow=n, ncol=n.bins)
  m.C[m.ind] <- p
  cs.C <- colSums(m.C)
  m.Y <- matrix(y, ncol=1)
  cs.Y <- colSums(m.Y)

  ## sets of indices for which colSums(A) > 0, resp. > 1
  ## in order to avoid division by zero
  d0 <- which(cs.A < 1)
  d1 <- which(cs.A < 2)
  #function that returns x_i : i \in d0
  d0.set <- function(x) {
    if (length(d0) > 0) {
      x <- x[-d0]
    }
    x
  }
  #function that returns x_i : i \in d1
  d1.set <- function(x) {
    if (length(d1) > 0) {
      x <- x[-d1]
    }
    x
  }

  ## calculate estimators
  rel <- 1/n * sum( d0.set((cs.B - cs.C)^2 / cs.A) )
  res <- 1/n * sum( d0.set(cs.A * (cs.B / cs.A - cs.Y / n)^2) )
  unc <- cs.Y * (n - cs.Y) / (n * n)

  ## correction factors for bias removal
  corr.s <- 1 / n * sum(d1.set(cs.B * (cs.A - cs.B) / (cs.A * (cs.A - 1))))
  corr.t <- cs.Y * (n - cs.Y) / n / n / (n-1)

  # avoid rel<0, res<0, res>1, and unc>1/4
  alpha <- min(c(rel/corr.s, 
                 max(c(res / (corr.s - corr.t), 
                       (res - 1) / (corr.s - corr.t))),
                 (1 - 4 * unc) / (4 * corr.t),
                 1)
              )

  rel2 <- rel - alpha * corr.s
  res2 <- res - alpha * corr.s + alpha * corr.t
  unc2 <- unc + alpha * corr.t
                  


  ###################################################
  # variances of reliability, resolution, uncertainty
  # estimated by propagation of uncertainty
  ###################################################

  m.cent <- diag(n) - 1 / n

  #####
  # REL
  #####
  d.rel.d.a <- -1 * (cs.B - cs.C) * (cs.B - cs.C) / (n * cs.A * cs.A)
  d.rel.d.b <- 2 * (cs.B - cs.C) / (n * cs.A)
  d.rel.d.c <- -2 * (cs.B - cs.C) / (n * cs.A)
  if (length(d0) > 0) {
    d.rel.d.a[d0] <- 0
    d.rel.d.b[d0] <- 0
    d.rel.d.c[d0] <- 0
  }
  jacobian.rel <- matrix(c(d.rel.d.a, d.rel.d.b, d.rel.d.c), nrow=1)

  m.X <- cbind(m.A, m.B, m.C)
  cov.X.rel <- t(m.X) %*% m.cent %*% m.X

  var.rel <- drop( jacobian.rel %*% cov.X.rel %*% t(jacobian.rel) )

  #####
  # RES
  #####
  d.res.d.a <- -1 / n * (cs.B / cs.A - cs.Y / n) * 
                      (cs.B / cs.A + cs.Y / n)
  d.res.d.b <- 2 / n * (cs.B / cs.A - cs.Y / n)
  d.res.d.y <- 0
  if (length(d0) > 0) {
    d.res.d.a[d0] <- 0
    d.res.d.b[d0] <- 0
  }
  jacobian.res <- matrix(c(d.res.d.a, d.res.d.b, d.res.d.y), nrow=1)

  m.X <- cbind(m.A, m.B, m.Y)
  cov.X.res <- t(m.X) %*% m.cent %*% m.X
  
  var.res <- drop( jacobian.res %*% cov.X.res %*% t(jacobian.res) )

  #####
  # UNC
  #####
  var.unc <- drop((1 - 2 * cs.Y / n) * (1 - 2 * cs.Y / n) / n / n * 
                   t(m.Y) %*% m.cent %*% m.Y)


  ######
  # REL'
  ######
  d.rel2.d.a <- -1 * ((cs.B - cs.C) * (cs.B - cs.C) + 
                      (cs.B * cs.B) / (cs.A - 1) - 
                      cs.A * cs.B * (cs.A - cs.B) / 
		      (cs.A - 1) / (cs.A - 1)
                     ) / (n * cs.A * cs.A)
  d.rel2.d.b <- (2 * cs.B - 1) / (n * cs.A - n) - 2 * cs.C / (n * cs.A)
  d.rel2.d.c <- -2 * (cs.B - cs.C) / (n * cs.A)
  if (length(d1) > 0) {
    d.rel2.d.a[d1] <- 0
    d.rel2.d.b[d1] <- 0
    d.rel2.d.c[d1] <- 0
  }
  jacobian.rel <- matrix(c(d.rel2.d.a, d.rel2.d.b, d.rel2.d.c), nrow=1)

  var.rel2 <- drop( jacobian.rel %*% cov.X.rel %*% t(jacobian.rel) )


  ######
  # RES'
  ######

  d.res2.d.a <- -1 / n * (cs.B / cs.A - cs.Y / n) * (cs.B / cs.A + cs.Y / n) +
                 cs.B / (n * cs.A * cs.A * (cs.A - 1) * (cs.A - 1)) *
                 ((cs.A - cs.B) * (cs.A - cs.B) - cs.B * (cs.B - 1))
  d.res2.d.b <- 2 / n * (cs.B / cs.A - cs.Y / n) - (cs.A - 2 * cs.B) /
                (n * cs.A * (cs.A - 1))
  d.res2.d.y <- (n - 2 * cs.Y) / n / n / (n - 1)
  if (length(d1) > 0) {
    d.res2.d.a[d1] <- 0
    d.res2.d.b[d1] <- 0
  }

  jacobian.res2 <- matrix(c(d.res2.d.a, d.res2.d.b, d.res2.d.y), nrow=1)

  var.res2 <- drop(jacobian.res2 %*% cov.X.res %*% t(jacobian.res2))

  ######
  # UNC'
  ######
 
  var.unc2 <- drop((n - 2 * cs.Y) * (n - 2 * cs.Y) / n / n / (n - 1) / (n - 1) *
                   t(m.Y) %*% m.cent %*% m.Y)
  


  #store everything
  bride.list <- list(
                  p = p,
                  y = y,
                  n.bins = n.bins,
                  rel = rel,
                  res = res,
                  unc = unc,
                  rel2 = rel2,
                  res2 = res2,
                  unc2 = unc2,
                  br = rel - res + unc,
                  rel.var = var.rel,
                  res.var = var.res,
                  unc.var = var.unc,
                  rel2.var = var.rel2,
                  res2.var = var.res2,
                  unc2.var = var.unc2
                )

  #return
  bride.list 
}

