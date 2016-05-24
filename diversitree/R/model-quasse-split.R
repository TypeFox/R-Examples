## 1: make
make.quasse.split <- function(tree, states, states.sd, lambda, mu,
                              nodes, split.t, control=NULL,
                              sampling.f=NULL) {
  cache <- make.cache.quasse.split(tree, states, states.sd,
                                   lambda, mu, nodes, split.t,
                                   control, sampling.f)
  n.part <- cache$n.part

  all.branches <- make.all.branches.quasse.split(cache, cache$control)
  rootfunc <- make.rootfunc.split(cache, make.rootfunc.quasse(cache))
  f.pars <- make.pars.quasse.split(cache)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.f=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.f, intermediates)
  }

  class(ll) <- c("quasse.split", "quasse", "dtlik", "function")
  ll
}

## 2: make.cache (& initial.tip)
make.cache.quasse.split <- function(tree, states, states.sd,
                                    lambda, mu, nodes, split.t,
                                    control, sampling.f) {
  cache <- make.cache.quasse(tree, states, states.sd, NULL, NULL,
                             control, sampling.f, TRUE)
  cache <- make.cache.split(tree, cache, nodes, split.t)

  n.part <- cache$n.part

  ## Speciation/extinction functions
  ## Most of the processing of these is currently done below.
  tmp.l <- check.f.quasse.split(lambda, n.part)
  tmp.m <- check.f.quasse.split(mu, n.part)
  i <- rbind(lambda=tmp.l$n, mu=tmp.m$n, drift=1, diffusion=1)
  j <- split(seq_len(sum(i)), rep(seq_along(i), i))
  dim(j) <- dim(i)
  rownames(j) <- c("lambda", "mu", "drift", "diffusion")

  names <- mapply(c,
                  lapply(tmp.l$names, sprintf, fmt="l.%s"),
                  lapply(tmp.m$names, sprintf, fmt="m.%s"),
                  "drift", "diffusion", SIMPLIFY=FALSE)
  argnames <- sprintf("%s.%d", unlist(names),
                      rep(seq_along(names), sapply(names, length)))

  cache$sampling.f <- check.sampling.f.split(sampling.f, 1, n.part)
  cache$aux.i <- seq_len(cache$control$nx)
  cache$lambda <- tmp.l$f
  cache$mu <- tmp.m$f
  cache$args <- t(j)
  
  cache$info$argnames <- argnames
  cache$info$lambda <- lambda
  cache$info$mu     <- mu

  cache
}

## This is actually a bad time sink; can take 0.01s from a 0.024s
## minimal evaluation.  Really, we don't need to compute this each
## time, but just subset by excluding nkl and nkr points.
initial.tip.quasse.split <- function(cache, control, x) {
  nx <- control$nx * control$r
  npad <- nx - length(x)
  tips <- cache$tips
  e0 <- lapply(cache$sampling.f, function(f)
               rep(1 - f, nx))[cache$group.nodes[tips]]
  y <- mapply(function(e0, mean, sd)
              c(e0, dnorm(x, mean, sd), rep(0, npad)),
              e0, cache$states, cache$states.sd, SIMPLIFY=FALSE)
  
  dt.tips.ordered(y, tips, cache$len[tips])
}

## This is almost identical to make.all.branches.split.dtlik except:
## 1. We make initial.conditions inside this function
## 2. No ode check
## 3. make.branches.quasse not make.branches.dtlik (also for aux)
## 4. The initial tip conditions are computed at the beginning of the
##    function.
make.all.branches.quasse.split <- function(cache, control) {
  control <- check.control.split(control)
  caching.branches <- control$caching.branches
  initial.conditions <- make.initial.conditions.quasse(control)

  branches.main <- make.branches.quasse(cache, control)
  branches.aux <- make.branches.aux.quasse(cache, control)
  branches.split <- make.branches.split(cache, branches.main,
                                        branches.aux, control)
  initial.conditions.split <-
    make.initial.conditions.split(cache, initial.conditions)

  function(pars, intermediates, preset=NULL) {
    if ( caching.branches )
      caching.branches.set.pars(pars, branches.split)
    cache$y <- initial.tip.quasse.split(cache, control, pars[[1]]$hi$x)
    all.branches.list(pars, cache,
                      initial.conditions.split,
                      branches.split, preset)
  }
}

make.branches.aux.quasse <- function(cache, control) {
  ## TODO: strictly, these should be backends...
  if ( control$method == "fftC" )
    branches.aux <- make.branches.aux.quasse.fftC(control, cache$sampling.f)
  else
    stop("Cannot do split QuaSSE with method ", control$method)
  branches.aux
}

make.pars.quasse.split <- function(cache) {
  n.args <- length(cache$info$argnames)
  n.part <- cache$n.part
  args <- cache$args
 
  function(pars) {
    names(pars) <- NULL # Because of use of do.call, strip names
    if ( length(pars) != n.args )
      stop(sprintf("Incorrect number of arguments (expected %d, got %d)",
                   n.args, length(pars)))
    drift <- pars[unlist(args[,3])]
    diffusion <- pars[unlist(args[,4])]

    ext <- quasse.extent(cache$control, drift, diffusion)

    ## expand the parameters, with our current extent.
    pars.l <- lapply(seq_len(n.part), function(i)
                     expand.pars.quasse(cache$lambda[[i]],
                                        cache$mu[[i]],
                                        args[i,], ext, pars))

    lambda.x <- unlist(lapply(pars.l, function(x) x$hi$lambda))
    mu.x <- unlist(lapply(pars.l, function(x) x$hi$mu))
    check.pars.quasse(lambda.x, mu.x, drift, diffusion)

    pars.l
  }
}

check.f.quasse.split <- function(f, rep) {
  if ( is.function(f) ) {
    n <- rep.int(check.f.quasse(f), rep)
    names <- rep.int(list(names(formals(f))[-1]), rep)
    f <- rep.int(list(f), rep)
  } else {
    if ( length(f) != rep )
      stop("Invalid length for speciation/extinction function")
    n <- sapply(f, check.f.quasse)
    names <- lapply(f, function(x) names(formals(x))[-1])
  }
  list(n=n, f=f, names=names)
}
