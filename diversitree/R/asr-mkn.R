make.asr.marginal.mkn <- function(lik, ...) {
  if ( inherits(lik, "mkn.ode") )
    stop("ASR not yet possible with ode-based Mkn")
  e <- environment(lik)
  f.pars <- e$f.pars
  all.branches <- e$all.branches
  cache <- e$cache
  
  info <- get.info(lik)
  k <- info$k
  states.idx <- info$idx.d
  
  cache.C <- list(parent=toC.int(cache$parent),
                  children=toC.int(t(cache$children)),
                  root=toC.int(cache$root))
  nodes.all.C <- toC.int(cache$root:max(cache$order))
  env <- new.env()

  function(pars, nodes=NULL, root=ROOT.FLAT, root.p=NULL, ...) {
    pars2 <- f.pars(pars)
    nodes.C <- if (is.null(nodes))
      nodes.all.C else toC.int(nodes + cache$n.tip)
    root.f <- function(pars, vals, lq)
      rootfunc.mkn(list(vals=vals, lq=lq), pars, root, root.p, FALSE)
    
    res <- all.branches(pars2, TRUE, NULL)
    .Call("r_asr_marginal_mkn", k, pars2, nodes.C, cache.C, res,
          root.f, env, PACKAGE="diversitree")
  }
}

make.asr.joint.mkn <- function(lik, ...) {
  k <- get.info(lik)$k
  cache <- get.cache(lik)
  is.mk2 <- inherits(lik, "mk2")

  order.C <- toC.int(rev(cache$order))
  parent.C <- toC.int(cache$parent)

  function(pars, n=1, ...) {
    ## TODO: this probably becomes something else...
    obj <- attr(lik(pars, intermediates=TRUE, ...), "intermediates")
    do.asr.joint(n, k, order.C, parent.C, obj$init, obj$pij,
                 obj$root.p, is.mk2)
  }
}

make.asr.stoch.mkn <- function(lik, slim=FALSE, ...) {
  is.mk2 <- inherits(lik, "mk2")
  cache <- get.cache(lik)
  k <- as.integer(get.info(lik)$k)

  if ( is.mk2 ) {
    states     <- as.integer(cache$states - 1L)
    states.pos <- c(0L, 1L)
  } else {
    states     <- as.integer(cache$states)
    states.pos <- seq_len(k)
  }

  edge.1 <- cache$edge[,1]
  edge.2 <- cache$edge[,2]
  edge.length <- cache$edge.length # == cache$len[edge.2]

  ## Joint distribution
  joint <- make.asr.joint(lik)

  ## Single branch simulator
  ptr <- .Call("r_smkn_alloc", k, 100L, PACKAGE="diversitree")

  function(pars, n=1, ...) {
    if ( n > 1 )
      stop("Not yet implemented (n>1)")

    ## 1: Draw a random sample from the nodes' joint distribution:
    node.state <- joint(pars, n, ...)

    ## If we are using mkn, we need to deflate the node states onto
    ## base-0 indices for the simulation code.
    if ( is.mk2 )
      anc.state <- as.integer(c(states, node.state))
    else
      anc.state <- as.integer(c(states - 1L, node.state - 1L))

    ## 2: Simulate branches.
    state.beg <- anc.state[edge.1]
    state.end <- anc.state[edge.2]
    history <- .Call("r_smkn_scm_run_all", ptr, pars, edge.length,
                     state.beg, state.end, is.mk2, slim,
                     PACKAGE="diversitree")

    if ( slim ) {
      ret <- list(node.state=node.state, history=history)
    } else {
      if ( !is.null(cache$node.label) )
        names(node.state) <- cache$node.label
      ret <- make.history(NULL, states, node.state, history,
                          TRUE, states.pos, check=FALSE)
    }
    ret
  }
}

## This appears totally broken at present, but it's not exported.
summarise.histories.mk2 <- function(x, phy) {
  summarise.branch <- function(x) {
    n <- length(x)
    j <- seq_len(n)
    y <- cbind(do.call(rbind, x),
               rep(j, sapply(x, nrow)), deparse.level=0)
    y <- y[order(y[,1]),]
    cbind(c(0, y[-j,1]),
          (c(0, cumsum(y[-j,2] * 2 - 1)) + sum(y[j,2])) / n,
          deparse.level=0)
  }

  nn <- length(x[[1]]$history)
  nx <- length(x)
  tmp <- matrix(unlist(lapply(x, "[[", "history"), FALSE), nn, nx)
  h <- apply(tmp, 1, summarise.branch)

  node.state <- rowMeans(matrix(unlist(lapply(x, "[[", "node.state")),
                                nn/2, nx))
  names(node.state) <- names(x[[1]]$node.state)

  states <- x[[1]]$tip.state
  make.history(phy, states, node.state, h, FALSE, 0:1)
}
