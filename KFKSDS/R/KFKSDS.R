
KFKSDS.deriv.C <- function(y, ss)
{
  n <- length(y)
  r <- ncol(ss$V)
  m <- ncol(ss$Z)
  nm <- n * m
  nr <- n * r
  rp1 <- r + 1
  nrp1 <- nr + n #n * rp1
  nrp1m <- nrp1 * m

  res <- .C("KFKSDS_deriv_C", 
    dim = as.integer(c(n, 1, ncol(ss$Z), r, 1, -1, 0)), 
    sy = as.numeric(y), sZ = as.numeric(ss$Z), sT = as.numeric(t(ss$T)), 
    sH = as.numeric(ss$H), sR = as.numeric(t(ss$R)), sV = as.numeric(ss$V), 
    sQ = as.numeric(ss$Q), sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
    dvof = double(nrp1), epshat = double(n), vareps = double(n),
    etahat = double(nr), vareta = double(nr),
    r = double(nm), N = double(nm), 
    dr = double(nrp1m), dN = double(nrp1m), 
    dahat = double(nrp1m), dvareps = double(nrp1), #dvarahat = double(nrp1m),
    PACKAGE = "KFKSDS")

  etahat <- matrix(res$etahat, nrow = n, ncol = r, byrow = TRUE)
  vareta <- matrix(res$vareta, nrow = n, ncol = r, byrow = TRUE)
  r <- matrix(res$r, nrow = n, ncol = m, byrow = TRUE)
  N <- matrix(res$N, nrow = n, ncol = m, byrow = TRUE)
  dr <- array(res$dr, dim = c(m, rp1, n))
  dN <- array(res$dN, dim = c(m, rp1, n))
  dahat <- array(res$dahat, dim = c(m, rp1, n))
  dvareps <- matrix(res$dvareps, nrow = n, ncol = rp1, byrow = TRUE)

  list(epshat = res$epshat, vareps = res$vareps,
    etahat = etahat, vareta = vareta,
    r = r, N = N, dr = dr, dN = dN, 
    dahat = dahat, dvareps = dvareps)
}

KFKSDS.deriv.steady.C <- function(y, ss, convergence = c(0.001, 10, 1.2))
{
  n <- length(y)
  r <- ncol(ss$V)
  m <- ncol(ss$Z)
  nm <- n * m
  nr <- n * r
  rp1 <- r + 1
  nrp1 <- nr + n #n * rp1
  nrp1m <- nrp1 * m

  res <- .C("KFKSDS_deriv_steady_C", 
    dim = as.integer(c(n, 1, ncol(ss$Z), r, 1, -1, 0)), 
    sy = as.numeric(y), sZ = as.numeric(ss$Z), sT = as.numeric(t(ss$T)), 
    sH = as.numeric(ss$H), sR = as.numeric(t(ss$R)), sV = as.numeric(ss$V), 
    sQ = as.numeric(ss$Q), sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
    tol = as.numeric(convergence[1]), maxiter = as.integer(convergence[2]), 
    ksconvfactor = as.numeric(convergence[3]),
    dvof = double(nrp1), epshat = double(n), vareps = double(n),
    etahat = double(nr), vareta = double(nr),
    r = double(nm), N = double(nm), 
    dr = double(nrp1m), dN = double(nrp1m), 
    dahat = double(nrp1m), dvareps = double(nrp1), #dvarahat = double(nrp1m),
    PACKAGE = "KFKSDS")

  etahat <- matrix(res$etahat, nrow = n, ncol = r, byrow = TRUE)
  vareta <- matrix(res$vareta, nrow = n, ncol = r, byrow = TRUE)
  r <- matrix(res$r, nrow = n, ncol = m, byrow = TRUE)
  N <- matrix(res$N, nrow = n, ncol = m, byrow = TRUE)
  dr <- array(res$dr, dim = c(m, rp1, n))
  dN <- array(res$dN, dim = c(m, rp1, n))
  dahat <- array(res$dahat, dim = c(m, rp1, n))
  dvareps <- matrix(res$dvareps, nrow = n, ncol = rp1, byrow = TRUE)

  list(epshat = res$epshat, vareps = res$vareps,
    etahat = etahat, vareta = vareta,
    r = r, N = N, dr = dr, dN = dN, 
    dahat = dahat, dvareps = dvareps)
}
