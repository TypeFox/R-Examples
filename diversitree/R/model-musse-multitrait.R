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

## TODO: I'm not sure if n.trait belongs in cache or in info?

## This differs from normal MuSSE in the interpretation of the
## arguments only; the actual calculation of the likelihood is done in
## the normal MuSSE functions.  However, there is also some changes to
## the initial conditions, etc, through the cache.
make.musse.multitrait <- function(tree, states, sampling.f=NULL,
                                  depth=NULL, allow.multistep=FALSE,
                                  strict=TRUE, control=list()) {
  cache <- make.cache.musse.multitrait(tree, states, depth,
                                       allow.multistep, sampling.f,
                                       strict)
  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions.musse)
  rootfunc <- rootfunc.musse
  f.pars.musse <- make.pars.musse.multitrait(cache)
  f.pars       <- make.pars.musse(cache$info$k)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE, pars.only=FALSE) {
    pars.musse <- f.pars.musse(pars)
    if ( pars.only )
      return(pars.musse)
    ## Below here is identical to as musse's likelihood function:
    pars2 <- f.pars(pars.musse)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("musse.multitrait", "musse", "dtlik", "function")
  ll
}

## 2:
make.info.musse.multitrait <- function(n.trait, depth, tr, phy) {
  k <- 2^n.trait
  ret <- make.info.musse(k, phy)
  ret$name <- "musse.multitrait"
  ret$name.pretty <- "MuSSE (multitrait)"
  ret$name.ode <- "musse"
  ret$mcmc.lowerzero <- FALSE
  ret$n.trait <- n.trait
  ret$argnames <- colnames(tr)
  ret
}

## 3: make.cache
##
## TODO: Handling sampling.f is not trivial; we'll assert that it is
## scalar for now.  It would be simple enough to allow the full 'k'
## long vector here, but people should not have to care about that.
##
## In theory if it was of length n.trait we could probably combine
## things appropriately, assuming multiplicative interactions.
## However, this is a weird assumption to have to make.
make.cache.musse.multitrait <- function(tree, states, depth=NULL,
                                        allow.multistep=FALSE,
                                        sampling.f=NULL,
                                        strict=TRUE) {
  tree <- check.tree(tree)
  states <- check.states.musse.multitrait(tree, states, strict=strict,
                                          strict.vals=0:1)
  n.trait <- ncol(states)
  k <- 2^n.trait
  depth <- check.depth(depth, n.trait)  

  cache <- make.cache(tree)
  cache$tr <- musse.multitrait.translate(n.trait, depth,
                                         colnames(states),
                                         allow.multistep)
  cache$info <- make.info.musse.multitrait(n.trait, depth, cache$tr,
                                           tree)
  cache$n.trait    <- n.trait
  cache$states     <- states
  cache$sampling.f <- rep(check.sampling.f(sampling.f, 1), k)
  cache$y <- initial.tip.musse.multitrait(cache)
  cache
}
initial.tip.musse.multitrait <- function(cache) {
  k <- cache$info$k
  n.trait <- cache$n.trait # == log2(k)
  f <- cache$sampling.f

  ## First, recode the states as a character code, replacing all NA
  ## values with '.' (which will match any state) and concatenating.
  code <- as.matrix(cache$states)
  storage.mode(code) <- "character"
  code[is.na(code)] <- "."
  code <- apply(code, 1, paste, collapse="")
  types <- unique(code)

  ## This is the "key"; the mapping from a series of binary traits
  ## onto the MuSSE mapping
  key <- apply(do.call(expand.grid, rep(list(0:1), n.trait)),
               1, paste, collapse="")

  y <- lapply(lapply(types, grep, key), function(i)
              c(1-f, ifelse(seq_len(k) %in% i, f, 0)))
  y.i <- match(code, types)

  dt.tips.grouped(y, y.i, cache)
}

######################################################################
## Additional functions
starting.point.musse.multitrait <- function(tree, lik, q.div=5,
                                            yule=FALSE) {
  if ( is.constrained(lik) ) {
    p <- starting.point.musse.multitrait(tree, attr(lik, "func"),
                                         q.div, yule)
    return(p[argnames(lik)])
  }

  if ( !inherits(lik, "musse.multitrait") )
    stop("'lik' must be a musse.multitrait model")
  
  tr <- get.cache(lik)$tr

  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  r <- if  ( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]

  p <- rep(0, ncol(tr))
  names(p) <- colnames(tr)
  p[c("lambda0", "mu0")] <- pars.bd
  p[grep("^q.+0$", names(p))] <- r / q.div
  p
}

## Parameter translation:
musse.multitrait.translate <- function(n.trait, depth=NULL,
                                       names=NULL,
                                       allow.multistep=FALSE) {
  if ( is.null(names) )
    names <- LETTERS[seq_len(n.trait)]
  depth <- check.depth(depth, n.trait)

  block.diagonal(reparam.1(n.trait, depth[1], names, "lambda"),
                 reparam.1(n.trait, depth[2], names, "mu"),
                 reparam.q(n.trait, depth[3], names, allow.multistep))
}

reparam.1 <- function(n.trait, depth, names, prefix=NULL) {
  X <- data.frame("0"=rep(1, length.out=2^n.trait), check.names=FALSE)
  if ( depth >= 1 ) {
    X <- cbind(X, do.call(expand.grid, rep(list(0:1), n.trait)))
    names(X)[-1] <- names

    if ( depth >= 2 ) {
      tmp <- unlist(lapply(2:depth, combn, x=names, simplify=FALSE),
                    FALSE)
      for ( i in tmp )
        X[[paste(i, collapse="")]] <- apply(X[i], 1, prod)
    }
    suffix.out <- apply(X[names], 1, paste, collapse="")
  } else {
    suffix.out <- apply(do.call(expand.grid, rep(list(0:1), n.trait)),
                        1, paste, collapse="")
  }

  X <- as.matrix(X)
  if ( !is.null(prefix) ) {
    colnames(X) <- paste(prefix, colnames(X), sep="")
    rownames(X) <- paste(prefix, suffix.out,  sep="")
  }

  X
}

reparam.q <- function(n.trait, depth, names, allow.multistep=FALSE) {
  k <- 2^n.trait
  from <- rep(seq_len(k), each=k-1)
  to <- unlist(lapply(seq_len(k), function(i) seq_len(k)[-i]))

  key <- do.call(expand.grid, rep(list(0:1), n.trait))
  names(key) <- names

  names.out <- sprintf("q%s.%s",
                       apply(key[from,,drop=FALSE], 1, paste, collapse=""),
                       apply(key[to,,drop=FALSE],   1, paste, collapse=""))

  nc <- rowSums(abs(key[from,,drop=FALSE] - key[to,,drop=FALSE]))
  n <- sum(choose(n.trait-1, 0:depth)) # (2^(n.trait-1))^(depth)
  X <- as.data.frame(matrix(0, nrow=length(from), ncol=2*n.trait*n))
  for ( i in seq_len(n.trait) ) {
    idx01 <- which(key[from,i] == 0 & key[to,i] == 1 & nc == 1)
    idx10 <- which(key[from,i] == 1 & key[to,i] == 0 & nc == 1)

    ## This is entirely incorrect, because n is incorrect...
    target <- seq((i - 1) * (2 * n) + 1, length.out=2 * n)
    tmp <- reparam.1(n.trait-1, depth, names[-i])
    X[idx01,target[1:n]] <- X[idx10,target[-(1:n)]] <- tmp
    names(X)[target] <- sprintf("q%s%s.%s", names[i],
                                rep(c("01", "10"), each=n),
                                colnames(tmp))
  }

  if ( allow.multistep ) {
    i <- nc > 1
    X2 <- matrix(0, length(nc), sum(i))
    X2[cbind(which(i), seq_len(sum(i)))] <- 1
    colnames(X2) <- names.out[i]
    X <- cbind(X, X2)
  }

  X <- as.matrix(X)
  rownames(X) <- names.out
  X
}

make.pars.musse.multitrait <- function(cache) {
  k <- cache$info$k
  n.trait <- cache$n.trait
  tr <- cache$tr
  npar <- ncol(tr)

  function(pars) {
    if ( length(pars) != npar )
      stop(sprintf("Invalid length parameters (expected %d)", 
                   npar))
    drop(tr %*% pars)
  }
}

############################################################
## Checking functions:
## Version of check.states that works with multitrait data properly.
check.states.musse.multitrait <- function(tree, states,
                                          allow.unnamed=FALSE,
                                          strict=FALSE,
                                          strict.vals=NULL) {
  if ( !is.data.frame(states) )
    stop("states must be a data.frame")
  if ( "0" %in% names(states) )
    stop("'0' cannot be in names (reserved for intercept)")
  if ( any(nchar(names(states)) != 1) )
    stop("All names must be length 1")
  if ( any(duplicated(names(states))) )
    stop("All names must be unique")
  
  if (is.null(rownames(states))) {
    if (allow.unnamed) {
      if (length(states) == length(tree$tip.label)) {
        rownames(states) <- tree$tip.label
        warning("Assuming states are in tree$tip.label order")
      }
      else {
        stop(sprintf("Invalid states length (expected %d)", 
                     length(tree$tip.label)))
      }
    }
    else {
      stop("The states matrix must contain rownames")
    }
  }

  if (!all(tree$tip.label %in% rownames(states))) 
    stop("Not all species have state information")

  if ( !is.null(strict.vals) ) {
    if ( isTRUE(all.equal(strict.vals, 0:1)) ) {
      i.log <- sapply(states, is.logical)
      states[,i.log] <- sapply(states[,i.log], as.integer)
    }
    
    if ( strict ) {
      f <- function(x)
        !isTRUE(all.equal(sort(strict.vals),
                          sort(unique(na.omit(x)))))
      if ( any(sapply(states, f)) )
        stop("Because strict state checking requested, all (and only) ", 
             sprintf("states in %s are allowed",
                     paste(strict.vals,  collapse=", ")))
    } else {
      f <- function(x)
        length(setdiff(sort(unique(na.omit(x))), strict.vals)) > 0
      if ( any(sapply(states, f)) )
        stop("Unknown states not allowed in states vector")
    }
  }

  as.matrix(states[tree$tip.label,,drop=FALSE])
}

check.depth <- function(depth, n.trait) {
  if ( is.null(depth) ) {
    depth <- c(n.trait, n.trait, n.trait-1)
  } else {
    depth <- check.par.length(depth, 3)
    depth[3] <- min(c(depth[3], n.trait-1))
  }
  if ( any(depth < 0) )
    stop("'depth' must be nonnegative")
  check.integer(depth)
  if ( any(depth > c(n.trait, n.trait, n.trait - 1)) )
    stop("requested depth too large")
  depth
}
