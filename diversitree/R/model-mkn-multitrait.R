## This differs from normal mkn in the interpretation of the arguments
## only; the actual calculation of the values is passed off to mkn.

## TODO: Updating this to allow missing data will require pulling it
## apart like make.musse.multitrait().  This should happen at some
## point.

## 1: make
make.mkn.multitrait <- function(tree, states, depth=NULL,
                                allow.multistep=FALSE,
                                strict=TRUE, control=list()) {
  cache <- make.cache.mkn.multitrait(tree, states, depth,
                                     allow.multistep, strict,
                                     control)
  lik.mkn <- make.mkn(tree, cache$states.mkn, cache$info$k,
                      strict=FALSE, control=cache$control.mkn)
  f.pars <- make.pars.mkn.multitrait(cache)

  ll <- function(pars, root=ROOT.OBS, root.p=NULL,
                 intermediates=FALSE, pars.only=FALSE) {
    pars2 <- f.pars(pars)
    if ( pars.only )
      return(pars2)
    lik.mkn(pars2, root, root.p, intermediates)
  }
  class(ll) <- c("mkn.multitrait", class(lik.mkn))
  ll
}

make.info.mkn.multitrait <- function(n.trait, tr, phy) {
  k <- 2^n.trait
  ret <- make.info.mkn(k, phy)
  ret$name <- "mkn.multitrait"
  ret$name.pretty <- "Mk(n) (multitrait)"
  ret$name.ode <- "mkn"
  ret$mcmc.lowerzero <- FALSE
  ret$n.trait <- n.trait
  ret$argnames <- colnames(tr)
  ret
}

make.cache.mkn.multitrait <- function(tree, states, depth,
                                      allow.multistep, strict,
                                      control) {
  n.trait <- ncol(states)
  k <- 2^n.trait
  
  if ( is.null(control$method) )
    control$method <- if (n.trait > 3) "ode" else "pij"

  states <- check.states.musse.multitrait(tree, states, strict=strict,
                                          strict.vals=0:1)
  ## TODO: This cannot be hard to do, if I can do it in MuSSE!
  ## The main issue is getting ode and exp to agree.
  if ( any(is.na(states)) )
    stop("Missing data not yet allowed")

  cache <- list()

  ## #### Make new states
  ## TODO/NEW:
  ## I think that this is actually somewhat general, and something
  ## very similar is used in initial.tip.musse.multitrait().
  ## I should make a file 'multitrait' and pull general support stuff
  ## like that into that file.
  code <- as.matrix(states)
  storage.mode(code) <- "character"
  code[is.na(code)] <- "."
  code <- apply(code, 1, paste, collapse="")
  types <- unique(code)

  ## This is the "key"; the mapping from a series of binary traits
  ## onto the Mk mapping
  key <- apply(do.call(expand.grid, rep(list(0:1), n.trait)),
               1, paste, collapse="")

  states.mkn <- match(code, key)
  names(states.mkn) <- rownames(states)
  cache$states.mkn <- states.mkn
  ## #### End

  cache$control.mkn <- control
  cache$tr <- mkn.multitrait.translate(n.trait, depth,
                                       colnames(states),
                                       allow.multistep)
  cache$info <- make.info.mkn.multitrait(n.trait, cache$tr, tree)  
  cache
}

mkn.multitrait.translate <- function(n.trait, depth=NULL,
                                     names=NULL,
                                     allow.multistep=FALSE) {
  if ( is.null(names) )
    names <- LETTERS[seq_len(n.trait)]

  if ( is.null(depth) )
    depth <- n.trait - 1
  else if ( length(depth) != 1 )
    stop("Depth must be 1")
  else if ( depth < 0 )
    stop("'depth' must be nonnegative")
  else if ( depth > n.trait - 1 )
    stop("requested depth too large")

  reparam.q(n.trait, depth, names, allow.multistep)
}

make.pars.mkn.multitrait <- function(cache) {
  n.trait <- cache$n.trait
  k <- cache$info$k
  tr <- cache$tr
  npar <- ncol(tr)
  f.pars <- make.pars.mkn(k)

  function(pars, pars.only=FALSE) {
    if ( length(pars) != npar )
      stop(sprintf("Invalid length parameters (expected %d)", 
                   npar))
    drop(tr %*% pars)
  }
}
