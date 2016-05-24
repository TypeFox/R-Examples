## Parameters come in the order:
##   l_1, ..., l_k, m_1, ..., m_k, qmat
## where qmat is the standard q matrix without the diagonals (same
## order as for make.musse).  As a result the parameters
## vector is 2k + k^2 - k = k(k+1) elements long.
make.tree.musse <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                            single.lineage=TRUE) {
  k <- (sqrt(1 + 4*length(pars))-1)/2
  if ( !isTRUE(all.equal(k, as.integer(k))) )
    stop("Invalid parameter length: must be k(k+1) long")
  check.pars.musse(pars, k)
  if ( x0 < 1 || x0 > k )
    stop("x0 must be an integer in [1,k]")

  pars <- cbind(matrix(pars[seq_len(2*k)], k, 2), 
                matrix(pars[(2*k+1):(k*(k+1))], k, k-1, TRUE))
  to <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])),
               k, k-1, TRUE)

  r.i <- rowSums(pars)

  extinct <- FALSE
  split   <- FALSE
  parent <- 0
  n.i <- rep(0, k)
  len <- 0
  t <- 0
  hist <- list()

  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0] <- 1
    start <- 0
  } else {
    ##states <- rep(x0, 2)
    ##n.taxa <- lineages <- n.i[x0] <- 2
    stop("Nope.")
  }

  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Pick a state, and then a lineage in that state, to modify, then
    ## an event type: 1: speciation, 2: extinction, >2: char change.
    state <- sample(k, 1, FALSE, r.n/r.tot)
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]
    type <- sample(k+1, 1, FALSE, pars[state,])

    if ( type == 1 ) { # Speciation:
      if ( n.taxa == max.taxa )
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      states[new.i] <- state
      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      n.i[state] <- n.i[state] + 1
      n.taxa <- n.taxa + 1
      lineages <- which(!split & !extinct)
    } else if ( type == 2 ) { # Extinction:
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1
    } else {
      states[lineage] <- state.new <- to[state,type - 2]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1,-1)
      hist[[length(hist)+1]] <- c(lineage, t, state, state.new)
    }
  }
  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     start=start, state=states, extinct=extinct,
                     split=split)

  hist <- as.data.frame(do.call(rbind, hist))
  if ( nrow(hist) == 0 )
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0
  
  attr(info, "t") <- t
  attr(info, "hist") <- hist
  info
}

tree.musse <- function(pars, max.taxa=Inf, max.t=Inf,
                       include.extinct=FALSE, x0=NA) {
  k <- (sqrt(1 + 4*length(pars))-1)/2
  if ( !isTRUE(all.equal(k, as.integer(k))) )
    stop("Invalid parameter length: must be k(k+1) long")
  if ( length(x0) != 1 || is.na(x0) || x0 < 1 || x0 > k )
    stop("Invalid root state")

  info <- make.tree.musse(pars, max.taxa, max.t, x0)
  ## This is a bit simple-minded, as extinction that erodes the root
  ## node will affect this once pruned.  However, it is always correct
  ## for trees without extinct taxa pruned.
  ## attr(info, "t") <- attr(info, "t") - info$len[1]
  phy <- me.to.ape.bisse(info[-1,], info$state[1])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

tree.musse.multitrait <- function(pars, n.trait, depth, max.taxa=Inf,
                                  max.t=Inf, include.extinct=FALSE,
                                  x0=NA) {
  tr <- musse.multitrait.translate(n.trait, depth)
  np <- ncol(tr)
  if ( !(all(is.finite(x0)) && length(x0) == n.trait) )
    stop(sprintf("x0 must be a vector of length %d"), n.trait)
  if ( !all(x0 %in% c(0, 1)) )
    stop("All x0 must be in [0,1]")
  if ( length(pars) != ncol(tr) )
    stop(sprintf("Expected %d parameters:\n\t%s",
                 ncol(tr), paste(colnames(tr), collapse=", ")))
  pars2 <- drop(tr %*% pars)

  code <- paste(x0, collapse="")
  types <- do.call(expand.grid, rep(list(0:1), n.trait))
  key <- apply(types, 1, paste, collapse="")
  x02 <- match(code, key)

  tree <- tree.musse(pars2, max.taxa, max.t, include.extinct, x02)
  
  tree$tip.state <- types[tree$tip.state,,drop=FALSE]
  rownames(tree$tip.state) <- tree$tip.label
  colnames(tree$tip.state) <- LETTERS[seq_len(n.trait)]
  tree
}
