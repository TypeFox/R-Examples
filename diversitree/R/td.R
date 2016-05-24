## Time chunks.

make.all.branches.td.dtlik <- function(cache, control,
                                       initial.conditions) {
  control <- check.control.ode(control)

  branches.td <- make.branches.td.dtlik(cache$info, control)
  initial.conditions.td <-
    make.initial.conditions.td(initial.conditions)
  function(pars, intermediates, preset=NULL)
    all.branches.matrix(pars, cache, initial.conditions.td,
                        branches.td, preset)
}

make.branches.td.dtlik <- function(info, control)
  make.branches.td(make.branches.dtlik(info, control))

make.branches.td <- function(branches) {
  function(y, len, pars, t0, idx) {
    t <- pars[,1]
    bpars <- pars[,-1]

    t1 <- t0 + len
    lq.carry <- 0

    i <- which(t > t0)[1]
    j <- which(t > t1[length(t1)])[1]

    done <- rep(FALSE, length(len))
    out <- vector("list", length(t))

    base <- matrix(NA, length(y), length(len))
    lq <- numeric(length(len))

    for ( epoch in i:j ) {
      not.last <- epoch < j
      ## Identify branches that terminate in this epoch:
      k <- t1 < t[epoch] & !done
      
      times <- c(t1[k], if (not.last) t[epoch])
      n <- length(times)

      res <- branches(y, times-t0, bpars[epoch,], t0, idx)
      res[[1]] <- res[[1]] + lq.carry

      if ( any(k) ) {
        base[,k] <- if ( not.last ) res[[2]][,-n] else res[[2]]
        lq[k]    <- if ( not.last ) res[[1]][-n]  else res[[1]]
        done[k] <- TRUE
      }
      
      if ( not.last ) {
        y <- res[[2]][,n]
        lq.carry <- res[[1]][n]
        t0 <- t[epoch]
      }
    }

    list(lq, base)
  }
}

make.initial.conditions.td <- function(initial.conditions) {
  function(init, pars, t, idx)
    initial.conditions(init, get.par.td(pars, t), t, idx)
}

make.rootfunc.td <- function(cache, rootfunc) {
  t.root <- cache$depth[cache$root]
  function(ans, pars, condition.surv, root, root.p, intermediates) {
    if ( root == ROOT.EQUI )
      stop(paste("Can't use ROOT.EQUI:", 
                 "stationary frequencies are messy for time models."))
    rootfunc(ans, get.par.td(pars, t.root), condition.surv, root,
             root.p, intermediates)
  }
}

## TODO/NEW: try this:
## make.rootfunc.t <- function(cache, rootfunc) {
##   t.root <- cache$depth[cache$root]
##   function(ans, pars, ...)
##     rootfunc(ans, get.par.td(pars, t.root), ...)
## }


update.info.td <- function(info, n.epoch) {
  n.epoch <- check.n.epoch(n.epoch)
  argnames.td <- c(sprintf("t.%d", seq_len(n.epoch-1)),
                   argnames.twopart(info$argnames, n.epoch))

  info$time.chunks <- TRUE
  info$argnames <- argnames.td
  info$name.ode <- info$name
  info$name.pretty <- sprintf("%s (time-chunks)", info$name.pretty)
  info$name <- sprintf("%s.td", info$name)  

  info$n.epoch <- n.epoch
  info
}

######################################################################
## Parameter generator function:
get.par.td <- function(pars, t)
  pars[which(pars[,1] > t)[1],-1]

######################################################################
## Checking
check.n.epoch <- function(n.epoch) {
  n.epoch <- check.integer(n.epoch)
  if ( !is.finite(n.epoch) || n.epoch < 1 )
    stop("n.epoch must be finite and at least 1")
  n.epoch
}

