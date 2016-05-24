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

## 1: make
make.bisse <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, strict=TRUE, control=list()) {
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  unresolved <- cache$unresolved
  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions.bisse)
  rootfunc <- rootfunc.musse

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    check.pars.bisse(pars)
    preset <- branches.unresolved.bisse(pars, unresolved)
    ans <- all.branches(pars, intermediates, preset)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse", "dtlik", "function")
  ll
}

## 2: info
make.info.bisse <- function(phy) {
  list(name="bisse",        # canonical name (as in make.<name>)
       name.pretty="BiSSE", # pretty name for printing
       ## Parameters:
       np=6L, # number of parameters the ODE takes, not the function.
       argnames=default.argnames.bisse(),
       ## Variables:
       ny=4L,     # number of variables
       k=2L,      # number of states (NA if continuous)
       idx.e=1:2, # index of 'E' variables, for xxsse models
       idx.d=3:4, # index of 'D' variables, for xxsse and mk models
       ## R version of derivative function:
       derivs=derivs.bisse,
       ## Phylogeny:
       phy=phy,   # here to help with printing, possibly plotting
       ## Inference:
       ml.default="subplex", # default ML search
       mcmc.lowerzero=TRUE,  # all paramters positive? (for mcmc)
       ## These are optional
       doc=NULL,    # extra string to print during print()
       reference=c( # vector of references
         "Maddison et al. (2007) doi:10.1080/10635150701607033",
         "FitzJohn et al. (2009) doi:10.1093/sysbio/syp067"))
}
default.argnames.bisse <- function() 
  c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")

## 3: make.cache (& initial.tip)
make.cache.bisse <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10,
                             strict=TRUE) {
  tree <- check.tree(tree)
  states <- check.states(tree, states,
                         strict=strict && is.null(unresolved),
                         strict.vals=0:1)

  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- make.unresolved.bisse(tree$clades, states)
    states <- states[tree$tip.label]
  }
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")

  cache <- make.cache(tree)
  cache$info <- make.info.bisse(tree)
  cache$states     <- states
  cache$sampling.f <- check.sampling.f(sampling.f, 2)

  cache$unresolved <- check.unresolved(cache, unresolved, nt.extra)
  if ( !is.null(cache$unresolved) ) {
    cache$tips   <- cache$tips[-cache$unresolved$i]
    cache$states <- cache$states[-cache$unresolved$i]
  }

  if ( strict && !is.null(unresolved) ) {
    seen <- ((0:1) %in% unique(na.omit(cache$states)) ||
             sapply(cache$unresolved[c("n0", "n1")], sum) > 0)
    if ( !all(seen) )
      stop("Because strict state checking requested, all (and only) ",
           "states in 0 1, are allowed")
  }

  if ( length(cache$tips) > 0 ) # in case all tips are unresolved.
    cache$y <- initial.tip.xxsse(cache, base.zero=TRUE)

  cache
}

## 4: initial.conditions:
## Note that we ignore both 't' and 'idx'.
initial.conditions.bisse <- function(init, pars, t, idx)
  c(init[c(1,2),1],
    init[c(3,4),1] * init[c(3,4),2] * pars[c(1,2)])

branches.unresolved.bisse <- function(pars, unresolved) {
  if ( is.null(unresolved) )
    return(NULL)
  Nc <- unresolved$Nc
  k <- unresolved$k
  nsc <- unresolved$nsc
  t <- unresolved$len
  nt <- max(Nc) + unresolved$nt.extra
  
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  base <- bucexpl(nt, lambda0, lambda1, mu0, mu1, q01, q10, t,
                  Nc, nsc, k)[,c(3,4,1,2),drop=FALSE]

  q <- rowSums(base[,3:4,drop=FALSE])
  base[,3:4] <- base[,3:4] / q

  ## Note the transpose here.
  list(target=unresolved$target,
       lq=log(q),
       base=t(base))
}

###########################################################################
## Additional functions
stationary.freq.bisse <- function(pars) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  g <- (lambda0 - mu0) - (lambda1 - mu1)
  eps <- (lambda0 + mu0 + lambda1 + mu1)*1e-14
  if ( abs(g) < eps ) {
    if ( q01 + q10 == 0 )
      p <- 0.5
    else
      p <- q10/(q01 + q10)
  } else {
    roots <- quadratic.roots(g, q10+q01-g, -q10)
    roots <- roots[roots >= 0 & roots <= 1]
    if ( length(roots) > 1 )
      p <- NA
    else
      p <- roots
  }
  c(p, 1-p)
}

starting.point.bisse <- function(tree, q.div=5, yule=FALSE) {
  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  if  ( pars.bd[1] > pars.bd[2] )
    p <- rep(c(pars.bd, (pars.bd[1] - pars.bd[2]) / q.div), each=2)
  else
    p <- rep(c(pars.bd, pars.bd[1] / q.div), each=2)
  names(p) <- default.argnames.bisse()
  p
}

check.unresolved <- function(cache, unresolved, nt.extra) {
  if ( is.null(unresolved) ) {
    return(NULL)
  } else if ( nrow(unresolved) == 0 ) {
    warning("Ignoring empty 'unresolved' argument")
    return(NULL)
  }

  required <- c("tip.label", "Nc", "n0", "n1")
  if ( !all(required %in% names(unresolved)) )
    stop("Required columns missing from unresolved clades")

  unresolved$tip.label <- as.character(unresolved$tip.label)
  if ( !all(unresolved$tip.label %in% cache$tip.label) )
    stop("Unknown tip species in 'unresolved'")

  if ( max(unresolved$Nc + nt.extra) > 200 )
    stop("The largest unresolved clade supported has %d species",
         200 - nt.extra)

  unresolved$i <- match(unresolved$tip.label, cache$tip.label)
  unresolved$target <- cache$tips[unresolved$i]

  unresolved$k   <- unresolved$n1
  unresolved$nsc <- unresolved$n0 + unresolved$n1
  unresolved <- as.list(unresolved)
  unresolved$nt.extra <- nt.extra

  unresolved$len <- cache$len[unresolved$target]

  unresolved
}

## For historical and debugging purposes, not used directly in the
## calculations, but branches function is generated this way
## internally.
make.branches.bisse <- function(cache, control)
  make.branches.dtlik(cache$info, control)

## This is here for reference, but not exported yet.  It should be
## tweaked in several ways
##   1. Starting parameter guessing should be done internally, at
##      least as an option.
##   2. Better listing of arguments
##   3. Automatic parsing of results into some sort of table; this
##      proabably requires classing this.
all.models.bisse <- function(f, p, ...) {
  f3 <- constrain(f, lambda1 ~ lambda0, mu1 ~ mu0, q01 ~ q10)
  f4.lm <- constrain(f, lambda1 ~ lambda0, mu1 ~ mu0)
  f4.lq <- constrain(f, lambda1 ~ lambda0, q01 ~ q10)
  f4.mq <- constrain(f, mu1 ~ mu0, q01 ~ q10)
  f5.l <- constrain(f, lambda1 ~ lambda0)
  f5.m <- constrain(f, mu1 ~ mu0)
  f5.q <- constrain(f, q01 ~ q10)

  ## Fit six and three parameter models
  if ( length(p) != 3 )
    stop("Starting point must be of length 3 (lambda, mu, q)")
  ans3 <- find.mle(f3, p, ...)

  ## Using the values from the 3p model, fit the 4p and 5p models:
  l <- ans3$par[1]
  m <- ans3$par[2]
  q <- ans3$par[3]

  ## Start the searches from the best model of the previous type.
  ans4.lm <- find.mle(f4.lm, c(l, m, q, q), ...)
  ans4.lq <- find.mle(f4.lq, c(l, m, m, q), ...)
  ans4.mq <- find.mle(f4.mq, c(l, l, m, q), ...)

  p.l <- if ( ans4.lm$lnLik > ans4.lq$lnLik )
    ans4.lm$par[c(1:2,2:4)] else ans4.lq$par[c(1:4,4)]
  p.m <- if ( ans4.lm$lnLik > ans4.mq$lnLik )
    ans4.lm$par[c(1,1:4)] else ans4.mq$par[c(1:4,4)]
  p.q <- if ( ans4.lq$lnLik > ans4.mq$lnLik )
    ans4.lq$par[c(1,1:4)] else ans4.mq$par[c(1:3,3:4)]
  ans5.l  <- find.mle(f5.l, p.l, ...)
  ans5.m  <- find.mle(f5.m, p.m, ...)
  ans5.q  <- find.mle(f5.q, p.q, ...)

  tmp <- list(ans5.l, ans5.m, ans5.q)
  i <- which.max(sapply(tmp, "[[", "lnLik"))
  p6 <- tmp[[i]]$par
  j <- list(c(1, 1:5), c(1:2, 2:5), c(1:5, 5))
  ans6 <- find.mle(f, p6[j[[i]]], ...)

  list(ans6=ans6,
       ans5.l =ans5.l,  ans5.m =ans5.m,  ans5.q =ans5.q,
       ans4.lm=ans4.lm, ans4.lq=ans4.lq, ans4.mq=ans4.mq,
       ans3=ans3)
}

check.pars.bisse <- function(pars)
  check.pars.nonnegative(pars, 6)

derivs.bisse <- function(t, y, pars) {
  E0 <- y[1]
  E1 <- y[2]
  D0 <- y[3]
  D1 <- y[4]

  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  c(-(mu0 + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1,
    -(mu1 + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0,
    -(mu0 + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1,
    -(mu1 + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0)
}
