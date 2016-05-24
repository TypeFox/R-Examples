## Split BD model.  This is basically MEDUSA, but a more concrete
## likelihood-function based version, rather than the highly optimised
## search version that Alfaro et al. describe.

## In the future, I hope to use the update() function to efficiently
## push this to match up with MEDUSA's capabilities.  As it stands,
## this will not be an efficient way of looping over nodes and running
## optimisations in the same way as MEDUSA.  But it should provide a
## decent reference implementation of the calculations.

## TODO: This function here only works with the "nee" method, and does
## not take a control argument.  It's also a complete mess.
make.bd.split <- function(tree, nodes, split.t=Inf, sampling.f=NULL,
                          unresolved=NULL) {
  cache <- make.cache.bd.split(tree, nodes, split.t, sampling.f,
                               unresolved)
  n.part <- cache$n.part

  all.branches <- make.all.branches.bd.split(cache)
  rootfunc <- make.rootfunc.bd.split(cache)

  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    pars <- check.pars.bd.split(pars, n.part)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, intermediates)
  }
  class(ll) <- c("bd.split", "bd", "dtlik", "function")
  ll
}

make.info.bd.split <- function(phy, nodes)
  update.info.split(make.info.bd(phy), nodes)

make.cache.bd.split <- function(tree, nodes, split.t=Inf,
                                sampling.f=NULL, unresolved=NULL) {
  tree <- check.tree(tree, node.labels=TRUE)

  if ( !isTRUE(all.equal(unique(split.t), Inf, check.attributes=FALSE)) )
    stop("split.t cannot yet be changed")
  nodes <- check.split(tree, nodes, split.t)$nodes
  n <- length(nodes)

  ## Process unresolved:
  unresolved <- check.unresolved.bd(tree, unresolved)

  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, n)

  ## Unresolved here is a different format to that expected by bisse,
  ## so this differs a little:
  edge <- tree$edge
  n.tip <- length(tree$tip.label)
  bt <- as.numeric(branching.times(tree))

  i <- seq_len(max(edge))
  j <- match(i, edge[,2])
  z <- cbind(anc=edge[j,1], dec=i,
             t.0=NA,
             t.1=bt[match(edge[j,1], (n.tip+1):max(edge))],
             t.len=tree$edge.length[j],
             n0=1, nt=NA,
             group=NA)
  z[n.tip + 1,1] <- n.tip + 1 # Special index for root.
  z[,"t.0"] <- z[,"t.1"] - z[,"t.len"]

  if ( is.null(unresolved) ) {
    z[match(seq_len(n.tip), z[,2]),"nt"] <- 1    
  } else {
    n.taxa <- unresolved$n[match(tree$tip.label, names(unresolved$n))]
    n.taxa[is.na(n.taxa)] <- 1
    z[match(seq_len(n.tip), z[,2]),"nt"] <- n.taxa
  }

  split.info <- make.cache.split(tree, list(), nodes[-1], split.t)
  z[,"group"] <- split.info$group.branches

  obj <- list(z=z, n.taxa=n.tip, n.node=tree$Nnode,
              sampling.f=sampling.f, t.root=max(bt),
              g.root=z[n.tip + 1,"group"], n.part=n)
  obj$info <- make.info.bd.split(tree, nodes)
  obj
}

make.all.branches.bd.split <- function(cache) {
  n.part <- cache$n.part
  ll.part <- lapply(seq_len(n.part), make.bd.split.part, cache=cache)

  function(pars, intermediates, preset=NULL) {
    res <- numeric(n.part)
    for ( i in seq_len(n.part) )
      res[i] <- ll.part[[i]](pars[[i]])
    list(vals=sum(res))
  }
}

make.rootfunc.bd.split <- function(cache) {
  ## I think that I can get this a lot easier somehow.
  z <- cache$z
  f <- cache$sampling.f[cache$g.root]
  t.root <- cache$t.root
  n <- tabulate(z[!is.na(z[,"nt"]),"group"], cache$n.part)
  root.constant <- 
    lfactorial(cache$n.taxa - 1) + sum(n*log(cache$sampling.f))

  function(vals, pars, condition.surv, intermediates) {
    if ( intermediates )
      stop("Sorry -- can't produce intermediates")
    loglik <- vals[[1]]

    if ( condition.surv ) {
      pars.r <- pars[[cache$g.root]]
      lambda <- pars.r[[1]]
      mu <- pars.r[[2]]
      r <- lambda - mu
      a <-  mu/lambda
      loglik <- loglik - log(f * f * r * (1 - a)) +
          2*log(abs(exp(-r * t.root)*(a - 1 + f) - f))
    }
    loglik + root.constant
  }
}

make.bd.split.part <- function(cache, i) {
  z <- cache$z[cache$z[,"group"] == i,]
  n.node <- sum(is.na(z[,"nt"]))

  z <- z[!is.na(z[,"t.1"]),] # drop root node

  f <- cache$sampling.f[i]

  t0 <- z[,"t.0"]
  t1 <- z[,"t.1"]
  dt <- z[,"t.len"]

  unresolved <- z[which(z[,"nt"] > 1),,drop=FALSE]
  
  function(pars) {
    lambda <- pars[1]
    mu <- pars[2]
    r <- lambda - mu

    ## The abs() here is because log(x^2) -> 2 log(abs(x))
    d <- r * dt +
      2*(log(abs((f * exp(r * t0) + 1-f) * lambda - mu)) -
         log(abs((f * exp(r * t1) + 1-f) * lambda - mu)))

    loglik <- sum(d) + n.node * log(lambda)

    if ( nrow(unresolved) > 0 ) {
      a <- mu / lambda
      ert <- exp(r * unresolved[,"t.1"])
      loglik <- loglik +
        sum((unresolved[,"nt"]-1) * (log(abs(ert - 1)) - log(abs(ert - a))))
    }

    loglik
  }
}

check.pars.bd.split <- function(pars, n.part) 
  check.pars.multipart(check.nonnegative(pars), n.part, 2)
