check.control.continuous <- function(control) {
  defaults <- list(method="vcv", backend="R")
  control <- modifyList(defaults, control)

  if ( length(control$method) != 1 )
    stop("control$method must be a scalar")
  if ( identical(control$method, "direct") ) {
    warning('method="direct" is deprecated, please use method="pruning"')
    control$method <- "pruning"
  }

  ## NOTE: contrasts is not actually working for everything yet (eb
  ## and ou models), and specifying it here will cause the function to
  ## silently select the pruning method.  But it's not the cutting
  ## edge if nobody is bleeding.
  methods <- c("vcv", "pruning", "contrasts")
  if ( !(control$method %in% methods) )
    stop(sprintf("control$method must be in %s",
                 paste(methods, collapse=", ")))

  if ( length(control$backend) != 1 )
    stop("control$backend must be a scalar")
  backends <- c("C", "R")
  if ( !(control$backend %in% backends) )
    stop(sprintf("control$backend must be in %s",
                 paste(backends, collapse=", ")))
  
  control
}

make.all.branches.continuous <- function(cache, control) {
  cache.C <- toC.cache.continuous(cache)

  br.name <- sprintf("branches_%s", cache$info$name)
  br <- getNativeSymbolInfo(br.name)$address
  ## Hard coded now, as bm and ou share.
  ic <- getNativeSymbolInfo("initial_conditions_bm")$address
  
  ptr <- .Call("r_make_dt_obj_cont", cache.C, ic, br,
               PACKAGE="diversitree")
  
  function(pars, intermediates=FALSE, preset=NULL) {
    if ( !is.null(preset) )
      stop("Don't know how to deal with preset values yet")
    res <- .Call("r_all_branches_cont", ptr, pars,
                 PACKAGE="diversitree")
    names(res) <- c("lq", "vals")
    if ( intermediates ) {
      vals <- res$vals
      res <- .Call("r_get_vals_cont", ptr, PACKAGE="diversitree")
      names(res) <- c("init", "base", "lq")
      res$vals <- vals
    }
    res    
  }
}

toC.cache.continuous <- function(cache) {
  ## This is super awkward for now, but may change in future.
  ## Basically, for this version, I don't want the indices ordered at
  ## all, and assert that the order should be 1..ntips.
  if ( length(cache$y$target) != cache$n.tip )
    stop("Missing tips")
  i <- order(cache$y$target)
  cache$y$target <- toC.int(cache$y$target[i]) # now seq_len(n.tip)-1L
  cache$y$y      <- cache$y$y[,i,drop=FALSE]
  cache$y$t      <- cache$y$t[i]
  cache$y$type   <- NULL

  cache$children <- toC.int(t(cache$children))
  cache$parent   <- toC.int(cache$parent)
  cache$order    <- toC.int(cache$order)
  cache$root     <- toC.int(cache$root)
  cache$n.tip    <- as.integer(cache$n.tip)
  cache$tips     <- toC.int(cache$tips)

  cache
}

make.all.branches.rescale.vcv <- function(cache, control) {
  n.tip <- cache$n.tip
  states <- cache$states
  states.sd <- cache$states.sd

  rescale <- make.rescale.phylo(cache$info$phy, cache$info$name)

  function(pars, intermediates) {
    s2 <- pars[[1]]

    vcv <- vcv.phylo(rescale(pars[-1]))
    vv <- s2 * vcv
    diag(vv) <- diag(vv) + states.sd^2

    VI <- solve(vv)
    mu <- sum(colSums(VI) * states) / sum(VI)
    dmvnorm2(states, rep(mu, n.tip), vv, VI, log=TRUE)
  }
}

make.all.branches.rescale.contrasts <- function(cache, control) {
  if (any(cache$states.sd > 0))
    stop("Cannot (yet) do contrasts based bm models with state error")

  ## This is annoying; we'll have to do this in make.bm.  If we don't,
  ## then we will reorder every time and that's going to hurt,
  ## timewise.  If I rewrite the pic() calculations to use my
  ## ordering, we can skip this step.  Given that eventually I think I
  ## want access to more of the pic components, that seems like a good
  ## idea.
  tree <- reorder(cache$info$phy, "pruningwise")
  states <- cache$states[tree$tip.label] # does reorder change this?
  rescale <- make.rescale.phylo(tree, cache$info$name)
  n <- length(tree$tip.label)

  function(pars, intermediates, preset=NULL) {
    s2 <- pars[[1]]
    tree <- rescale(pars[-1])

    ## Copied from model-bm.R:make.cache.bm
    ## There is duplication here with make.cache.pgls; perhaps merge
    ## these?  That might help with some of the root treatment
    ## things.
    pics <- pic(cache$states, tree, var.contrasts=TRUE)
    u <- pics[,"contrasts"]
    V <- pics[,"variance"]
    V0 <- pgls.root.var.bm(tree)
    ## This step is brutal.  We could rewrite the branch lengths there
    ## once the lookups are done properly.  This is particularly bad
    ## because it will involve doing another tree reordering.
    ## Costly!  What's more is we can actually get it from the PIC
    ## calculation.  So this part needs to be done separately, and for
    ## all the pic methods.  This entire section of code is horribly
    ## duplicated!
    root.x <- pgls.root.mean.bm(tree, cache$states)

    ## This is the bit that is shared with all.branches.bm
    ll <- -(n * log(2 * pi * s2) +
            sum(log(V)) +
            log(V0) +
            sum(u * u) / s2) * 0.5
    list(loglik=ll,
         root.x=root.x,
         root.v=s2 * V0,
         # not sure if these are needed...
         V=V, V0=V0)
  }
}
