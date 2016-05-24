make.all.branches.split.dtlik <- function(cache, control,
                                          initial.conditions) {
  control <- check.control.ode(control)
  control <- check.control.split(control)
  caching.branches <- control$caching.branches

  branches.main <- make.branches.dtlik(cache$info, control)
  branches.aux <- make.branches.aux.dtlik(cache, control)
  branches.split <- make.branches.split(cache, branches.main,
                                        branches.aux, control)
  initial.conditions.split <-
    make.initial.conditions.split(cache, initial.conditions)
  function(pars, intermediates, preset=NULL) {
    if ( caching.branches )
      caching.branches.set.pars(pars, branches.split)
    all.branches.matrix(pars, cache,
                        initial.conditions.split,
                        branches.split, preset)
  }    
}

## Convert a branches function into one that automatically accounts
## for splits.
## TODO: Deal with the case of a null branches.aux function so that
## this will work for bd-ode and mkn-ode.
make.branches.split <- function(cache, branches, branches.aux,
                                control) {
  group <- cache$group.branches
  split.with.edge <- cache$split.with.edge
  aux.use <- cache$aux.use
  aux.i <- cache$aux.i
  vars.in.matrix <- !is.na(cache$info$ny)

  if ( control$caching.branches ) {
    branches.c     <- make.caching.branches(cache, branches)
    branches.aux.c <- make.caching.branches.aux(cache, branches.aux)
  } else {
    branches.c <- branches
    branches.aux.c <- branches.aux
  }

  function(y, len, pars, t0, idx) {
    g <- group[idx]

    if ( length(len) == 1 ) {
      ## cat(sprintf("%d: %d, len=%2.5f, t0=%2.5f\n", idx, g, len, t0))
      ## 1: This might be an internal edge, so it is possible that
      ## there are auxilliary variables that need computing.
      ## (this could also be a terminal edge where there is only a
      ## single group)
      p <- pars[[g]]

      if ( !aux.use[idx] ) {
        ## Normal case: just run the branch
        branches.c(y, len, p, t0, idx)
      } else if ( split.with.edge[g] ) {
        ## This is a join.  In this one here, the edge also has the
        ## derived character state:
        ans <- branches.c(y, len, p, t0, idx)

        g.parent <- group[cache$parent[idx]]
        aux <- branches.aux.c(g.parent, t0+len, pars[[g.parent]], idx)
        if ( is.list(ans) )
          ans[[2]][aux.i] <- aux
        else
          ans[-1][aux.i] <- aux

        ans # return
      } else {
        ## Also a join, but now the edge is in the ancestral character
        ## state:
        y[aux.i] <- branches.aux.c(g, t0, p, idx)
        branches.c(y, len, p, t0, idx)
      }
    } else if ( length(unique(g)) == 1 ) {
      ## 2: This must be a set of tips, but it is only of one type, so
      ## let's do it simply.
      branches.c(y, len, pars[[g[1]]], t0, idx)
    } else {
      ## 3: This is a set of tips that have more than one group type
      ## in it, so loop over the different group types.
      grp <- sort.default(unique.default(group[idx]))
      i <- split.default(seq_along(idx), group[idx])

      if ( vars.in.matrix ) {
        lq <- numeric(length(len))
        ans <- matrix(NA, cache$info$ny, length(len))

        for ( g.i in seq_along(grp) ) {
          j <- i[[g.i]]
          tmp <- branches.c(y, len[i[[g.i]]], pars[[grp[g.i]]], t0, idx[j])
          lq[j] <- tmp[[1]]
          ans[,j] <- tmp[[2]]
        }
        list(lq, ans)
      } else {
        ans <- vector("list", length(len))
        for ( g.i in seq_along(grp) ) {
          j <- i[[g]]
          ans[j] <- branches.c(y, len[i[[g]]], pars[[grp[g]]], t0, idx[j])
        }
        ans
      }
    }
  }
}

## Compute E values (or other auxiliary variables) at a split.
## Unlike the other make.branches functions, this one needs the full
## cache so that we can get sampling.f out.  This is a bit of a pity.
make.branches.aux.dtlik <- function(cache, control) {
  ## The check will happen twice here, but that's OK.  It makes
  ## name.ode, which is useful.
  info.aux <- check.info.ode(cache$info, control)
  info.aux$name.ode <- sprintf("%s_aux", info.aux$name.ode)
  info.aux$ny       <- info.aux$k
  info.aux$idx.d    <- integer(0)
  info.aux$derivs   <- make.derivs.aux(info.aux$derivs, cache$aux.i, info.aux$k)
  branches <- make.branches.dtlik(info.aux, control)

  y <- lapply(cache$sampling.f, function(x) 1-x)
  n <- length(y)

  function(i, len, pars, idx) {
    if ( i > n )
      stop("No such partition")
    branches(y[[i]], len, pars, 0)[[2]]    
  }
}

make.initial.conditions.split <- function(cache, initial.conditions) {
  group <- cache$group.nodes
  function(init, pars, t, idx)
    initial.conditions(init, pars[[group[idx]]], t, idx)
}

## TODO/NEW: See the time-based models for what needs doing here to
## allow non-SSE models to have splits.
make.rootfunc.split <- function(cache, rootfunc) {
  grp.root <- cache$group.nodes[cache$root] # should always be 1?
  function(ans, pars, condition.surv, root, root.p, intermediates) {
    if ( root == ROOT.EQUI )
      stop(paste("Can't use ROOT.EQUI:", 
                 "stationary frequencies are messy for split models."))
    rootfunc(ans, pars[[grp.root]], condition.surv, root,
             root.p, intermediates)
  }
}

## TODO/NEW: try replacing with this:
## make.rootfunc.split <- function(cache, rootfunc) {
##   grp.root <- cache$group.nodes[cache$root] # should always be 1?
##   function(ans, pars, ...)
##     rootfunc(ans, pars[[grp.root]], ...)    
## }

## Modify a cache for use with a split model.
## 
## Determine the group membership of every edge.
## Augment the usual cache vector with some extra information:
make.cache.split <- function(phy, cache, nodes, split.t) {
  tmp <- split.group(phy, nodes, split.t)
  nodes <- tmp$nodes

  cache$nodes <- tmp$nodes
  cache$group.nodes <- tmp$group.nodes
  cache$group.branches <- tmp$group.branches
  cache$split.with.edge <- tmp$split.with.edge

  cache$n.part <- length(nodes)
  cache$aux.use <- rep(FALSE, length(cache$group.branches))
  cache$aux.use[nodes[-1]] <- TRUE

  cache$info <- update.info.split(cache$info, nodes)

  cache
}

update.info.split <- function(info, nodes) {
  n.part <- length(nodes)

  ## To catch an issue with bd.split
  if (is.null(info$name)) {
    info$name <- info$name.pretty
  }
  info$partitioned <- TRUE
  info$argnames <- argnames.twopart(info$argnames, n.part)
  info$name.ode <- info$name
  info$name.pretty <- sprintf("%s (split tree)", info$name.pretty)
  info$name <- sprintf("%s.split", info$name)

  info$n.part <- n.part
  info$nodes <- nodes

  info
}

update.info.uneven <- function(info, info.single) {
  info$name <- sprintf("%s.uneven", info.single$name)
  info$name.pretty <- sprintf("%s (uneven sampling tree)",
                              info.single$name.pretty)
  info$argnames <- info.single$argnames
  info
}

## As above, but do additional work with initial conditions and
## sampling.f for the xxSSE models.
make.cache.split.xxsse <- function(phy, cache, nodes, split.t,
                                   sampling.f) {
  cache <- make.cache.split(phy, cache, nodes, split.t)

  info <- cache$info
  n.part <- cache$n.part
  
  cache$sampling.f <-
    check.sampling.f.split(sampling.f, info$k, n.part)
  cache$aux.i <- info$idx.e
  cache$y <- make.initial.tip.xxsse.split(cache)

  cache
}

######################################################################
## Utility
split.group <- function(phy, nodes, split.t) {
  tmp <- check.split(phy, nodes, split.t)
  nodes <- tmp$nodes
  split.with.edge <- tmp$split.with.edge
  n.tip <- length(phy$tip.label)
  edge <- phy$edge

  descendants.idx <- function(node)
    which(edge[,1] == node | edge[,2] %in% descendants(node, edge))
  
  f <- function(node, phy, group) {
    if ( is.na(node) )
      return(list(group=group, base=NA))
    i <- which(edge[,1] == node | edge[,2] %in% descendants(node, edge))
    base <- group[edge[,2] == node]
    group[i[group[i] == base]] <- max(group) + 1
    list(group=group, base=base)
  }

  group <- rep(1, nrow(edge))

  base <- integer(length(nodes))
  for ( i in seq_along(nodes) ) {
    tmp <- f(nodes[i], phy, group)
    group <- tmp$group
    base[i] <- tmp$base
  }

  group <- group[match(seq_len(max(edge)), edge[,2])]
  group[n.tip + 1] <- 1

  group.branches <- group
  i <- which(!split.with.edge)
  if ( length(i) > 0 )
    group.branches[which(i)] <- base[which(i)]

  list(nodes=c(n.tip+1, nodes[-1]),
       group.nodes=group, group.branches=group.branches,
       split.with.edge=split.with.edge)
}

make.initial.tip.xxsse.split <- function(cache) {
  k <- cache$info$k
  n.part <- cache$n.part
  grp <- cache$group.branches  
  idx.D <- (k+1):(2*k)

  out <- vector("list", length(cache$y))
  for ( i in seq_along(cache$y)) {
    y.i <- cache$y[[i]]
    out[[i]] <- vector("list", n.part)
    
    for ( group in seq_len(n.part) ) {
      keep <- grp[y.i$target] == group
      if ( any(keep) ) {
        s.f <- cache$sampling.f[[group]]
        out[[i]][[group]] <- 
          list(y     = c(1-s.f, s.f * y.i$y[idx.D]),
               y.i   = y.i$y.i,
               target= y.i$target[keep],
               t     = y.i$t[keep])
      }
    }
  }
  unlist(out, FALSE)
}

######################################################################
## Caching branches

## The first of these that I'll do is fully passive.  It checks to
## see if the parameters and input are identical.  This should be
## the most basic, and most foolproof.  It may not be the fastest
## though, as some of these checks add overhead.  The use of
## identical should keep things fairly fast though.

## I want to make caching versions of the branches function.
caching.branches.set.pars <- function(p, branches, branches.aux) {
  e <- environment(branches)
  e.b <- environment(e$branches.c)
  e.a <- environment(e$branches.aux.c)

  e.a$pars.same <- e.b$pars.same <- i <-
    mapply(identical, p, e.b$prev.pars)
  e.a$prev.pars <- e.b$prev.pars <- p

  i
}

make.caching.branches <- function(cache, branches) {
  cat("Using experimental caching branches!\n")
  
  n <- length(cache$len)
  n.part <- cache$n.part
  group <- cache$group.branches
  vars.in.matrix <- !is.na(cache$info$ny)  

  pars.same <- logical(n.part)
  prev.pars <- vector("list", n.part)
  prev.lq <- numeric(n)
  
  if ( vars.in.matrix ) {
    prev.vars <- prev.base <- matrix(NA, cache$info$ny, n)
    branches.caching <- function(y, len, p, t0, idx) {
      g <- group[idx[1]]
      if ( pars.same[g] && identical(prev.vars[,idx[1]], y) ) {
        list(prev.lq[idx],
             prev.base[,idx,drop=FALSE])
      } else {
        prev.vars[,idx] <<- y
        ret <- branches(y, len, p, t0, idx)
        prev.lq[idx]    <<- ret[[1]]
        prev.base[,idx] <<- ret[[2]]
        ret
      }
    }
  } else { # Just QuaSSE, really.
    prev.vars <- prev.base <- vector("list", n)
    branches.caching <- function(y, len, p, t0, idx) {
      g <- group[idx[1]]
      if ( length(idx) > 1 ) # Will take some work, but not used atm.
        stop("vector idx not yet allowed")
      if ( pars.same[g] && identical(prev.vars[[idx]], y) ) {
        c(prev.lq[idx],
          prev.base[[idx]])
      } else {
        prev.vars[[idx]] <<- y
        ret <- branches(y, len, p, t0, idx)
        prev.lq[idx]   <<- ret[1]
        prev.base[[idx]] <<- ret[-1]
        ret
      }
    }
  }

  branches.caching
}

## Note that for now, the returned function here includes an
## additional argument to the input function:
##   normal aux: (g, t, p)
##   returned:   (g, t, p, idx)
make.caching.branches.aux <- function(cache, branches.aux) {
  n <- length(cache$len)
  prev.aux <- vector("list", n)
  pars.same <- logical(cache$n.part)

  branches.aux.caching <- function(g, t, p, idx) {
    if ( pars.same[g] ) {
      prev.aux[[idx]]
    } else {
      prev.aux[[idx]] <<- ret <- branches.aux(g, t, p)
      ret
    }
  }

  branches.aux.caching
}

make.derivs.aux <- function(derivs, idx, k) {
  force(derivs)
  y.extra <- rep(0, k)
  function(t, y, pars)
    derivs(t, c(y, y.extra), pars)[idx]
}

######################################################################
## Checking
check.split <- function(phy, nodes, split.t) {
  if ( length(split.t) == 1 && length(nodes) > 1 ) {
    if ( split.t == 0 || split.t == Inf )
      split.t <- rep(split.t, length(nodes))
    else
      stop("If split.t is length 1, it must be '0' or 'Inf'")
  } else if ( length(split.t) != length(nodes) ) {
    stop("'nodes' and 'split.t' must be the same length")
  }
    
  ## Check that all nodes are ok
  n.tip <- length(phy$tip.label)
  if ( is.character(nodes) )
    nodes <- match(nodes, phy$node.label) + n.tip
  if ( any(is.na(nodes) | nodes <= n.tip+1 | nodes > n.tip + phy$Nnode) )
    stop("Invalid node specification")
  
  ## Check that split times are OK
  partial <- !(split.t == Inf | split.t == 0)
  if ( any(partial) )
    stop("Partial split times not yet allowed")

  list(nodes=c(NA, nodes),
       split.with.edge=c(NA, split.t == Inf))
}

check.control.split <- function(control) {
  if ( is.null(control$caching.branches) )
    control$caching.branches <- FALSE
  else {
    val <- control$caching.branches
    if ( !(length(val) == 1 && val %in% c(TRUE, FALSE)) )
      stop("Invalid value for control$caching.branches")
  }
  control
}
