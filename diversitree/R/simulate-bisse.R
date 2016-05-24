## This is written in an odd style because I want to shift this into C
## soon after getting it working.
## TODO: I think the single.lineage bit is broken, so I am not
## exposing it in the next function up (the two-lineage case is
## broken).  This will interact particularly with the
##   me.to.ape.bisse(info[-1])
## line in tree.bisse() above.
make.tree.bisse <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                            single.lineage=TRUE) {
  pars <- matrix(pars, 2, 3)

  extinct <- FALSE
  split   <- FALSE
  parent <- 0
  n.i <- c(0, 0)
  r.i <- rowSums(pars)
  len <- 0
  t <- 0
  hist <- list()
  
  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0+1] <- 1
    start <- 0    
  } else {
    ##states <- rep(x0, 2)
    ##n.taxa <- lineages <- n.i[x0+1] <- 2
    stop("Nope.")    
  }

  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    ## When does an event happen?
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

    ## Proceed.  What state does an event happen to?
    state <- as.integer(runif(1) > r.n[1]/r.tot)
    state.i <- state + 1

    ## Pick a lineage for that state:
    j <- sample(n.i[state.i], 1)
    lineage <- lineages[states[lineages] == state][j]

    type <- sample(3, 1, FALSE, pars[state.i,])

    if ( type == 1 ) {
      ## Speciating:
      if ( n.taxa == max.taxa )
        ## Don't add this one.
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      states[new.i] <- state
      parent[new.i] <- lineage
      start[new.i] <- t      
      len[new.i] <- 0

      n.i[state.i] <- n.i[state.i] + 1
      n.taxa <- n.taxa + 1

      lineages <- which(!split & !extinct)

    } else if ( type == 2 ) {
      ## Extinct
      extinct[lineage] <- TRUE

      lineages <- which(!split & !extinct)

      n.i[state.i] <- n.i[state.i] - 1
      n.taxa <- n.taxa - 1
    } else {
      ## Character switch:
      n.i <- n.i + if ( state == 0 ) c(-1,1) else c(1,-1)
      states[lineage] <- 1 - state
      hist[[length(hist)+1]] <- c(lineage, t, state, 1-state)
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

tree.bisse <- function(pars, max.taxa=Inf, max.t=Inf,
                       include.extinct=FALSE, x0=NA) {
  ## This checking is duplicated in make.tree.bisse.C.core, for some
  ## reason.
  check.pars.bisse(pars)
  if ( is.na(x0) )
    x0 <- as.integer(runif(1) > stationary.freq.bisse(pars)[1])
  else if ( length(x0) != 1 || !(x0 == 0 || x0 == 1) )
    stop("Invalid root state")
  
  info <- make.tree.bisse(pars, max.taxa, max.t, x0)
  phy <- me.to.ape.bisse(info[-1,], info$state[1])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

make.tree.bisse.C <- function(pars, max.taxa, max.t, x0, n=100) {
  stop("I belive that this is broken and should not be used.")
  ans <- make.tree.bisse.C.core(pars, max.taxa, max.t, x0, n)

  i <- seq_len(ans$n.entries)
  info <- list(idx=i,
               len=ans$len[i],
               parent=ans$parent[i] + 1,
               start=ans$start[i],
               state=ans$states[i],
               extinct=as.logical(ans$extinct[i]),
               split=as.logical(ans$split[i]))
  attr(info, "rownames") <- i
  class(info) <- "data.frame"

  hist.t <- ans$hist.t[ans$hist.t > 0]
  if ( length(hist.t) > 0 ) {
    hist <- cbind(t(matrix(ans$hist[seq_len(length(hist.t) * 3)], 3)),
                  hist.t, deparse.level=0)[,c(1,4,2,3),drop=FALSE]
    hist <- as.data.frame(hist)
    names(hist) <- c("idx", "t", "from", "to")
    hist$idx <- hist$idx+1
    hist$x0 <- info$start[match(hist$idx, info$idx)]
    hist$tc <- hist$t - hist$x0
  } else {
    hist <- NULL
  }

  attr(info, "t") <- ans$t.start
  attr(info, "hist") <- hist  
  info
}

## TODO: There are a number of ways that will throw the C code here
## into an infinite loop; for instance, zero lambdas and mus (with
## nonzero qs) and max.t=Inf, etc.
## TODO: Still need to be able to handle max.t=Inf, before changing
## this for the main simulation code.
make.tree.bisse.C.core <- function(pars, max.taxa, max.t, x0, n=100,
                                   verbose=FALSE) {
  check.pars.bisse(pars)
  check.scalar(x0)
  if ( !is.finite(x0) || x0 < 0 || x0 > 1 )
    stop("x0 must be in [0,1] and non-NA")
  if ( is.infinite(max.taxa) && is.infinite(max.t) )
    stop("at most one of max.taxa and max.t may be infinite")
       
  n <- as.integer(n)
  pars2 <- t(matrix(pars, 2, 3))
  max.taxa <- as.integer(max.taxa)

  parent <- states <- extinct <- split <- integer(n)
  hist <- integer(3*n)
  start <- len <- hist.t <- numeric(n)
  n.entries <- as.integer(1)
  n.entries.hist <- as.integer(0)
  n.info <- as.integer(c(n.entries, n.entries.hist))

  parent[1] <- as.integer(-1)
  states[1] <- as.integer(x0)
  t <- 0.0

  verbose <- as.integer(verbose)

  ## Because this is an exponential growth, 10 iterations gets very
  ## large!  Starting at 100 and doubling each time 
  for ( attempt in 1:10 ) {
    ans <- .C("r_simulate_bisse", pars=pars2, max.taxa=max.taxa,
              max.t=max.t, parent=parent, states=states,
              extinct=extinct, split=split, start=start, len=len,
              hist=hist, hist.t=hist.t, n.info=n.info,
              n=n, t=t, verbose=verbose)
    if ( ans$t > 0 )
      break

    n <- as.integer(n * 2)
    i <- seq_len(n)
    j <- seq_len(n*3)

    parent <- states <- extinct <- split <- integer(n)
    hist <- integer(3*n)
    start <- len <- hist.t <- numeric(n)

    parent[i] <- ans$parent
    states[i] <- ans$states
    extinct[i] <- ans$extinct
    split[i] <- ans$split
    start[i] <- ans$start
    len[i] <- ans$len
    hist[j] <- ans$hist
    hist.t[i] <- ans$hist.t
    n.info <- ans$n.info
    t <- -ans$t
  }

  if ( ans$t < 0 )
    stop("This somehow failed...")

  i <- seq_len(ans$n.info[1])
  n.h <- ans$n.info[2]
  info <- list(idx=i,
               len=ans$len[i],
               parent=ans$parent[i] + 1,
               start=ans$start[i],
               state=ans$states[i],
               extinct=as.logical(ans$extinct[i]),
               split=as.logical(ans$split[i]),
               hist=ans$hist[seq_len(n.h*3)],
               hist.t=ans$hist.t[seq_len(n.h)],
               t=ans$t)
}

make.tree.bisse.R.core <- function(pars, max.taxa=Inf, max.t=Inf,
                                   x0) {
  pars <- matrix(pars, 2, 3)

  extinct <- FALSE
  split   <- FALSE
  parent <- 0
  n.i <- c(0, 0)
  r.i <- rowSums(pars)
  len <- 0
  t <- 0
  hist <- list()
  
  states <- x0
  n.taxa <- lineages <- n.i[x0+1] <- 1
  start <- 0    

  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    ## When does an event happen?
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

    ## Proceed.  What state does an event happen to?
    state <- as.integer(runif(1) > r.n[1]/r.tot)
    state.i <- state + 1

    ## Pick a lineage for that state:
    j <- sample(n.i[state.i], 1)
    lineage <- lineages[states[lineages] == state][j]

    type <- sample(3, 1, FALSE, pars[state.i,])

    if ( type == 1 ) {
      ## Speciating:
      if ( n.taxa == max.taxa )
        ## Don't add this one.
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      states[new.i] <- state
      parent[new.i] <- lineage
      start[new.i] <- t      
      len[new.i] <- 0

      n.i[state.i] <- n.i[state.i] + 1
      n.taxa <- n.taxa + 1

      lineages <- which(!split & !extinct)

    } else if ( type == 2 ) {
      ## Extinct
      extinct[lineage] <- TRUE

      lineages <- which(!split & !extinct)

      n.i[state.i] <- n.i[state.i] - 1
      n.taxa <- n.taxa - 1
    } else {
      ## Character switch:
      n.i <- n.i + if ( state == 0 ) c(-1,1) else c(1,-1)
      states[lineage] <- 1 - state
      hist[[length(hist)+1]] <- c(lineage, t, state, 1-state)
    }
  }

  info <- list(idx=seq_along(extinct),
               len=len,
               parent=parent,
               start=start,
               state=states,
               extinct=extinct,
               split=split,
               hist=hist,
               t=t)
}

make.tree.bisse.R <- function(pars, max.taxa, max.t, x0) {
  ans <- make.tree.bisse.R.core(pars, max.taxa, max.t, x0)

  info <- ans[1:7]
  attr(info, "row.names") <- info$idx
  class(info) <- "data.frame"

  hist <- as.data.frame(do.call(rbind, ans$hist))
  if ( nrow(hist) == 0 )
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0
  
  attr(info, "t") <- ans$t
  attr(info, "hist") <- hist  
  info
}
