## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches

# The full ClaSSE parameter structure is a vector with, for n states:
#     n * n * (n+1) / 2  speciation rates (lambda_ijk, j <= k)
#     n                  extinction rates (mu_i)
#     n * n - n          transition rates (q_ij, i != j)
#   = (n + 3) * n^2 / 2  elements

## 1: make
make.classe <- function(tree, states, k, sampling.f=NULL, strict=TRUE,
                        control=list()) {
  ## Note that this uses MuSSE's cache...
  cache <- make.cache.classe(tree, states, k, sampling.f, strict)
  initial.conditions <- make.initial.conditions.classe(k)
  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions)
  rootfunc <- rootfunc.classe
  f.pars <- make.pars.classe(k)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    ## TODO: This is different to other functions, as the
    ## stationary.freq function assumes the non-expanded case.
    ## However, it would be straightforward to modify stationary.freq
    ## classe to use the expanded case.  At worst, it could just strip
    ## off the extra parameters, but I think that it builds these
    ## anyway.
    ##
    ## This will be an issue for creating a split or time version.
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("classe", "dtlik", "function")
  ll
}

## 2: info
make.info.classe <- function(k, phy) {
  list(name="classe",
       name.pretty="ClaSSE",
       ## Parameters:
       np=as.integer((k+3)*k*k/2 + k),
       argnames=default.argnames.classe(k),
       ## Variables:
       ny=as.integer(2*k),
       k=as.integer(k),
       idx.e=as.integer(1:k),
       idx.d=as.integer((k+1):(2*k)),
       ## R version of derivatives function
       derivs=derivs.classe,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Goldberg (submitted)"))
}
default.argnames.classe <- function(k) {
  fmt <- sprintf("%%0%dd", ceiling(log10(k + .5)))
  sstr <- sprintf(fmt, 1:k)
  lambda.names <- sprintf("lambda%s%s%s", rep(sstr, each=k*(k+1)/2),
                          rep(rep(sstr, times=seq(k,1,-1)), k), 
                          unlist(lapply(1:k, function(i) sstr[i:k])))
  mu.names <- sprintf("mu%s", sstr)
  q.names <- sprintf("q%s%s", rep(sstr, each=k-1), 
                     unlist(lapply(1:k, function(i) sstr[-i])))
  c(lambda.names, mu.names, q.names)
}

## 3: make.cache (& initial.tip)
## Note: classe uses the same functions as musse here to generate the
## cache and initial tip states.
make.cache.classe <- function(tree, states, k, sampling.f=NULL,
                             strict=TRUE) {
  if (k > 31)
    stop("No more than 31 states allowed.  Increase in classe-eqs.c.")
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  cache$info <- make.info.classe(k, tree)
  cache
}

## 4: initial.conditions:
## save on index computations by wrappping initial.conditions.classe
make.initial.conditions.classe <- function(n) {
  ## n = number of states; called k elsewhere but k is used as an index below
  nseq <- seq_len(n)
  lam.idx <- matrix(seq_len(n*n*(n+1)/2), byrow=TRUE, nrow=n)

  idxD <- (n+1):(2*n)
  j <- rep(nseq, times=seq(n,1,-1))
  k <- unlist(lapply(1:n, function(i) nseq[i:n]))
  d <- rep(NA, n)

  initial.conditions.classe <- function(init, pars, t, is.root=FALSE) {
    ## E_i(t), same for N and M
    e <- init[nseq,1]

    ## D_i(t), formed from N and M
    DM <- init[idxD,1]
    DN <- init[idxD,2]
    DM.DN <- 0.5 * (DM[j] * DN[k] + DM[k] * DN[j])
    for (i in nseq)
      d[i] <- sum(pars[lam.idx[i,]] * DM.DN) # slower with apply

    ## a touch slower but cleaner:
    ##   idxlam = seq_len(n*n*(n+1)/2)
    ##     d = colSums(matrix(pars[idxlam], ncol=n) * DM.DN)
    ## or slightly better (but still not faster than for):
    ##   lamseq = seq_len(n*n*(n+1)/2)
    ##   lam.mat = matrix(lamseq, ncol=n)
    ##     lam.mat[lamseq] = pars[lamseq]
    ##     d = colSums(lam.mat * DM.DN)

    c(e, d)
  }
}

rootfunc.classe <- function(res, pars, condition.surv, root, root.p,
                            intermediates) {
  vals <- res$vals
  lq <- res$lq
  k <- length(vals)/2

  i <- seq_len(k)
  d.root <- vals[-i]

  ## TODO: This could be tidied up:
  root.equi <- function(pars) stationary.freq.classe(pars, k)
  root.p <- root.p.calc(d.root, pars, root, root.p, root.equi)
 
  if ( condition.surv ) {
    ## species in state i are subject to all lambda_ijk speciation rates
    nsum <- k*(k+1)/2
    lambda <- colSums(matrix(pars[1:(nsum*k)], nrow=nsum))
    e.root <- vals[i]
    d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2)
  }

  if ( root == ROOT.ALL )
    loglik <- log(d.root) + sum(lq)
  else
    loglik <- log(sum(root.p * d.root)) + sum(lq)

  if ( intermediates ) {
    res$root.p <- root.p
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- vals
  }

  loglik
}

###########################################################################
## Additional functions
## Heuristic starting point
## based on starting.point.geosse()
starting.point.classe <- function(tree, k, eps=0.5) {
  if (eps == 0) {
    s <- (log(Ntip(tree)) - log(2)) / max(branching.times(tree))
    x <- 0
    d <- s/10
  } else {
    n <- Ntip(tree)
    r <- ( log( (n/2) * (1 - eps*eps) + 2*eps + (1 - eps)/2 *
           sqrt( n * (n*eps*eps - 8*eps + 2*n*eps + n))) - log(2)
         ) / max(branching.times(tree))
    s <- r / (1 - eps)
    x <- s * eps
    q <- s - x
  }
  p <- c( rep(s / (k*(k+1)/2), k*k*(k+1)/2 ), rep(x, k), rep(q, k*(k-1)) )
  names(p) <- default.argnames.classe(k)
  p
}

stationary.freq.classe <- function(pars, k) {
  if (k == 2) {
    g <- (sum(pars[1:3]) - pars[7]) - (sum(pars[4:6]) - pars[8])
    eps <- sum(pars[1:8]) * 1e-14
    ss1 <- pars[9]  + 2*pars[3] + pars[2]  # shift from 1
    ss2 <- pars[10] + 2*pars[4] + pars[5]  # shift from 2

    if ( abs(g) < eps ) {
      if (ss1 + ss2 == 0) 
        eqfreq <- 0.5
      else
        eqfreq <- ss2/(ss1 + ss2)
      eqfreq <- c(eqfreq, 1 - eqfreq)
    } else {
      roots <- quadratic.roots(g, ss2 + ss1 - g, -ss2)
      eqfreq <- roots[roots >= 0 & roots <= 1]
      if ( length(eqfreq) > 1 )
        eqfreq <- NA
      else
        eqfreq <- c(eqfreq, 1 - eqfreq)
    }
  } else { ## also works for k=2, but much slower
    eqfreq <- stationary.freq.classe.ev(pars, k)
  }
  eqfreq
}

## like stationary.freq.geosse()
stationary.freq.classe.ev <- function(pars, k) {
  A <- projection.matrix.classe(pars, k)
  ## continuous time, so the dominant eigenvalue is the largest one
  ## return its eigenvector, normalized
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  evA$vectors[,i] / sum(evA$vectors[,i])
}

projection.matrix.classe <- function(pars, k) {
  A <- matrix(0, nrow=k, ncol=k)

  nsum <- k*(k+1)/2
  kseq <- seq_len(k)
  pars.lam <- pars[seq(1, nsum*k)]
  pars.mu <- pars[seq(nsum*k+1, (nsum+1)*k)]
  pars.q <- pars[seq((nsum+1)*k+1, length(pars))]

  ## array indices of lambda's in parameter vector
  idx.lam <- cbind(rep(kseq, each=nsum), rep(rep(kseq, times=seq(k,1,-1)), k),
                   unlist(lapply(kseq, function(i) i:k)))
  ## transpose of matrix indices of q's in parameter vector
  idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])), 
                 rep(kseq, each=k-1))

  ## take care of off-diagonal elements
  for (n in seq_len(nsum*k)) {
    ## add this lambda to A[daughter states, parent state]
    ## (separate steps in case the daughter states are the same)
    r <- idx.lam[n,]
    A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
    A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
  }
  A[idx.q] <- A[idx.q] + pars.q

  ## fix the diagonal elements
  diag(A) <- 0
  diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) 
                 sum(pars.lam[seq((i-1)*nsum+1, i*nsum)]) - pars.mu[i]))
  A
}

## For historical and debugging purposes, not used directly in the
## calculations, but branches function is generated this way
## internally.
make.branches.classe <- function(cache, control)
  make.branches.dtlik(cache$info, control)

## Parameter manipulation
make.pars.classe <- function(k) {
  np0 <- as.integer((k+3)*k*k/2)  # number of actual params
  np <- np0 + k                   # number, including Q's diagonal elements
  
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  x <- k * k * (k + 1) / 2 + k
  idx.lm <- seq_len(x)
  idx.q <- seq(x+1, np0)
  
  function(pars) {
    check.pars.classe(pars, k)
    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    c(pars[idx.lm], qmat)
  }
}

## These two functions are intended to make the classe parameters easier to
## visualize and populate, since they get unwieldy with more than two states.
## The speciation rate array is indexed lambda[parent state, daughter1 state,
## daughter2 state].  The transition matrix is indexed q[from state, to state].
## Elements that are not parameters get NA: daughter2 > daughter 1, from = to.
## The parameter list might be a good way to work with constrain(), eventually.

## Input: list containing lambda_ijk array, mu vector, q_ij array, num states
## Output: parameter vector, ordered as default.argnames.classe describes
flatten.pars.classe <- function(parlist) {
  k <- parlist$nstates
  kseq <- seq_len(k)

  idx.lam <- cbind( rep(kseq, each=k*(k+1)/2), 
                    rep(rep(kseq, times=seq(k,1,-1)), k), 
                    unlist(lapply(kseq, function(i) i:k)) )

  idx.q <- cbind( rep(kseq, each=k-1), 
                  unlist(lapply(kseq, function(i) (kseq)[-i])) )

  pars <- c(parlist$lambda[idx.lam], parlist$mu, parlist$q[idx.q])
  names(pars) <- default.argnames.classe(k)
  pars
}

## Output: list containing lambda_ijk array, mu vector, q_ij array, num states
## Input: parameter vector, ordered as default.argnames.classe describes
inflate.pars.classe <- function(pars, k) {
  check.pars.classe(pars, k)
  kseq <- seq_len(k)

  Lam <- array(NA, dim=rep(k, 3))  # 3 = parent + 2 daughters
  idx <- cbind(rep(kseq, each=k*(k+1)/2), rep(rep(kseq, times=seq(k,1,-1)), k),
               unlist(lapply(kseq, function(i) i:k)))
  j <- length(idx[,1])
  Lam[idx] <- pars[seq(j)]

  Mu <- pars[seq(j+1, j+k)]
  names(Mu) <- NULL

  Q <- array(NA, dim=rep(k, 2))
  idx <- cbind(rep(kseq, each=k-1), unlist(lapply(kseq, function(i) kseq[-i])))
  Q[idx] <- pars[seq(j+k+1, length(pars))]

  list(lambda=Lam, mu=Mu, q=Q, nstates=k)
}

check.pars.classe <- function(pars, k)
  check.pars.nonnegative(pars, (k+3)*k*k/2)

derivs.classe <- function(t, y, pars) {
  ## TODO: Need to write this (Emma: Do you have a copy somewhere?)
  stop("Not yet possible")
}
