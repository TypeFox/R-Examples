## Based on the contrasts approach in Freckleton 2012, MEE.
##
## This approach will not tolerate standard deviations for tip
## species.  Though incorporating them would probably not be all that
## hard.
##
## One way to look at that is how arbutus:::model.phylo.se() does this.

make.all.branches.bm.contrasts <- function(cache, control) {
  states <- cache$states
  if (any(cache$states.sd > 0))
    stop("Cannot (yet) do contrasts based bm models with state error")

  # Same names as make.all.branches.pgls.contrasts
  n      <- cache$n  # Number of tips
  u      <- cache$u  # Contrasts for the states
  V      <- cache$V  # Contrast variances
  V0     <- cache$V0 # Root contrast variance
  root.x <- cache$root.x

  ## Note that like all.branches.pgls there are constants here that
  ## can be factored out. n*log(2 * pi) + sum(log(V)) + log(V0) is
  ## constant.

  ## The sum(u * u) / sigma2 term is funny - because at ML, we take
  ## sum(u * u) / n -> sigma2, it turns out that sum(u*u)/s2 is n.

  ## The sweet thing about this approach is that it makes the
  ## connection clearer to what the different steps are doing.  If we
  ## can relate this to how the per-branch calculations happen in
  ## make.bm, then we're most of the way there.  My guess is that
  ## s2*log(V) will be partly what we're looking for here.
  
  function(pars, intermediates, preset=NULL) {
    s2 <- pars[[1]]
    ll <- -(n * log(2 * pi * s2) +
            sum(log(V)) +
            log(V0) +
            sum(u * u) / s2) * 0.5

    ## So, from the look of this; either I can just work with the way
    ## that I did this in PGLS, or I can return the list of three
    ## values.  To get this approximately correct, i'd have
    ##   vals[1] = root.x  -- stored
    ##   vals[2] = s2 * V0 -- computed
    ##   vals[3] = ll
    ## and then rootfunc.bm.pruning immediately works as currently
    ## implemented I think.  That also means that we can allow the
    ## root calculation to vary in OU type models.
    ##
    ## This is what I'm looking at modifying the root function by, but
    ## that won't quite give me what I want.
    ## dll <- -0.5 * (root.x.given - root.x)^2 / (s2 * V0)    
    ## as.numeric(dll + dll)
    list(loglik=ll,
         root.x=root.x,
         root.v=s2 * V0,
         # not sure if these are needed...
         V=V, V0=V0)
  }
}

## So, by default the calculations above seem to give us the
## appropriate answer for the case where the root is set to the ML
## value.  Which is nice.  But it means that some extra steps are
## required to apply different root treatments.
rootfunc.bm.contrasts <- function(res, pars, root, root.x,
                                  intermediates) {
  if (root == ROOT.MAX) {
    loglik <- res$loglik
  } else {
    root.v <- res$root.v
    loglik <- res$loglik + log(2 * pi * root.v) / 2
    if (root == ROOT.FLAT) {
      loglik <- loglik # pass
    } else if (root == ROOT.OBS) {
      loglik <- loglik - log(2 * sqrt(pi * root.v))
    } else if (root == ROOT.GIVEN) {
      if (is.null(root.x))
        stop("root.x not provided, but root=ROOT.GIVEN specified")
      loglik <- loglik + dnorm(root.x, res$root.x, sqrt(root.v), TRUE)
    } else {
      stop("Invalid root mode")
    }
  }

  if ( intermediates ) {
    res$root.p <- NA # not sure what would be good here...
    attr(loglik, "intermediates") <- res
  }

  loglik
}
