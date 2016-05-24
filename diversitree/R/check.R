## Checking utilities.  These are things that happen over and over
## again, and are tedious to have to write into each function.

## Only things that are not model specific should go in this file.  If
## the function ends in .xxx where 'xxx' is the name of a model, it
## probably belongs in model-xxx.R
check.tree <- function(tree, ultrametric=TRUE, bifurcating=TRUE,
                       node.labels=FALSE) {
  if ( !inherits(tree, "phylo") )
    stop("'tree' must be a valid phylo tree")
  if ( ultrametric && !is.ultrametric(tree) )
    stop("'tree' must be ultrametric")
  if ( any(tree$edge.length < 0) )
    stop("Negative branch lengths in tree")
  ## ape's is.binary.tree() can let a few nasties through - for
  ## e.g. each tritomy, an unbranched node and this gets through.
  ## This expression is a little stricter, even if a touch slower.
  if ( bifurcating && (!is.binary.tree(tree) ||
                       any(tabulate(tree$edge[, 1]) == 1)) )
    stop("'tree must be bifurcating (no polytomies or unbranched nodes)'")

  if ( any(duplicated(tree$tip.label)) )
    stop("Tree contains duplicated tip labels")
  
  if ( node.labels ) {
    if ( is.null(tree$node.label) )
      tree$node.label <- sprintf("nd%d", seq_len(tree$Nnode))
    else if ( any(duplicated(tree$node.label)) )
      stop("Tree contains duplicated node labels")
  }

  tree
}

check.states <- function(tree, states, allow.unnamed=FALSE,
                         strict=FALSE, strict.vals=NULL,
                         as.integer=TRUE) {
  if ( is.matrix(states) ) {
    ## Multistate characters (experimental).  This will not work with
    ## clade trees, but they are only interesting for BiSSE, which has
    ## NA values for multistate (even weight).
    if ( inherits(tree, "clade.tree") )
      stop("Clade trees won't work with multistate tips yet")
    n <- rowSums(states > 0)
    if ( any(n == 0) )
      stop(sprintf("No state found for taxa: %s",
                   paste(names(n)[n == 0], collapse=", ")))

    i.mono <- which(n == 1)
    i.mult <- which(n >  1)

    tmp <- matrix.to.list(states)
    names(tmp) <- rownames(states)

    states.mult <- lapply(tmp[i.mult], as.numeric)

    states <- rep(NA, length(tmp))
    names(states) <- names(tmp)
    states[i.mono] <- sapply(tmp[i.mono], function(x)
                             which(x != 0))

    attr(states, "multistate") <- list(i=i.mult, states=states.mult)
  }
  
  if ( is.null(names(states)) ) {
    if ( allow.unnamed ) {
      if ( length(states) == length(tree$tip.label) ) {
        names(states) <- tree$tip.label
        warning("Assuming states are in tree$tip.label order")
      } else {
        stop(sprintf("Invalid states length (expected %d)",
                     length(tree$tip.label)))
      }
    } else {
      stop("The states vector must contain names")
    }
  }
  
  if ( !all(tree$tip.label %in% names(states)) )
    stop("Not all species have state information")

  ## TODO: When multistate characters are present, this may fail even
  ## for cases where it should not.
  if ( !is.null(strict.vals) ) {
    if ( isTRUE(all.equal(strict.vals, 0:1)) )
      if ( is.logical(states) )
        states[] <- as.integer(states)
    
    if ( strict ) {
      if ( !isTRUE(all.equal(sort(strict.vals),
                             sort(unique(na.omit(states))))) )
        stop("Because strict state checking requested, all (and only) ",
             sprintf("states in %s are allowed",
                     paste(strict.vals, collapse=", ")))
    } else {
      extra <- setdiff(sort(unique(na.omit(states))), strict.vals)
      if ( length(extra) > 0 )
        stop(sprintf("Unknown states %s not allowed in states vector",
                     paste(extra, collapse=", ")))
    }
    if ( as.integer && any(!is.na(states)) )
      states <- check.integer(states)
  }

  if ( inherits(tree, "clade.tree") ) {
    spp.clades <- unlist(tree$clades)
    if ( !all(spp.clades %in% names(states)) )
      stop("Species in 'clades' do not have states information")
    states[union(tree$tip.label, spp.clades)]
  } else {
    ret <- states[tree$tip.label]
    ## Ugly hack...
    attr(ret, "multistate") <- attr(states, "multistate")
    ret
  }
}

check.par.length <- function(x, length) {
  if ( length(x) == 1 )
    rep(x, length)
  else if ( length(x) == length )
    x
  else
    stop(sprintf("'%s' of incorrect length",
                 deparse(substitute(x))))
}

check.sampling.f <- function(sampling.f, n) {
  if ( is.null(sampling.f) )
    sampling.f <- rep(1, n)
  else
    sampling.f <- check.par.length(sampling.f, n)

  if ( max(sampling.f) > 1 || min(sampling.f) <= 0 )
    stop("sampling.f must be on range (0,1]")
  sampling.f
}

check.sampling.f.split <- function(sampling.f, n, n.part) {
  if ( is.null(sampling.f) )
    rep(list(rep(1, n)), n.part)
  else if ( n == 1 && is.numeric(sampling.f) &&
           length(sampling.f == n.part) )
    check.sampling.f(sampling.f, n.part)
  else if ( is.numeric(sampling.f) )
    rep(list(check.sampling.f(sampling.f, n)), n.part)
  else if ( is.list(sampling.f) )
    lapply(sampling.f, check.sampling.f, n)
  else
    stop("Invalid sampling.f")
}

check.bounds <- function(lower, upper, x0=NULL) {
  if ( !is.null(x0) && (any(x0 < lower) || any(x0 > upper)) )
    stop("Starting parameter falls outside of problems bounds")
  if ( any(lower >= upper) )
    stop("'upper' must be strictly greater than 'lower'")
}

check.pars.multipart <- function(pars, n.part, n.per) {
  if ( is.matrix(pars) ) {
    if ( nrow(pars) != n.part )
      stop(sprintf("Expected %d parameter sets", n.part))
    if ( ncol(pars) != n.per )
      stop(sprintf("Expected %d parameters in each set", n.per))
    pars <- matrix.to.list(pars)
  } else if ( is.list(pars) ) {
    if ( length(pars) != n.part )
      stop(sprintf("Expected %d parameter sets", n.part))
    if ( !all(unlist(lapply(pars, length)) == n.per) )
      stop(sprintf("Expected %d parameters in each set", n.per))      
  } else {
    if ( length(pars) != n.part * n.per )
      stop(sprintf("Expected %d parameters", n.part * n.per))
    pars <- matrix.to.list(matrix(pars, n.part, n.per, TRUE))
  }
  pars
}

## Check that a number can reasonably be considered an integer.
check.integer <- function(x) {
  if ( is.null(x) )
    stop("NULL argument for ", deparse(substitute(x)))
  nna <- !is.na(x)
  if ( length(x) > 0 && !any(nna) )
    stop("No non-NA values for ", deparse(substitute(x)))
  if ( length(x) && max(abs(x[nna] - round(x[nna]))) > 1e-8 )
    stop("Non-integer argument for ", deparse(substitute(x)))
  storage.mode(x) <- "integer"
  x
}

check.scalar <- function(x) {
  if ( length(x) != 1 )
    stop(deparse(substitute(x)), " must be a scalar")
  x
}

check.control.ode <- function(control=list()) {
  defaults <- list(tol=1e-8, backend="gslode",
                   eps=0,                  # branch calculations
                   gsl.stepper="rkck",     # gslode specific
                   safe=TRUE, unsafe=FALSE # deSolve specific
                   )
  control <- modifyList(defaults, control)

  if ( !("compiled" %in% names(control)) )
    control$compiled <- control$backend == "gslode"

  backends <- c("deSolve", "gslode")
  if ( length(control$backend) != 1 )
    stop("'backend' must be a single element")
  control$backend <- backends[pmatch(control$backend, backends)]
  if ( is.na(control$backend) )
    stop("Invalid backend selected")

  if ( control$compiled && control$backend != "gslode" )
    stop("Compiled derivatives (control: compiled=TRUE) only with gslode")

  control$tol <- check.scalar(control$tol)
  control$eps <- check.scalar(control$eps)
  control$safe <- check.scalar(control$safe)

  if ( !is.numeric(control$tol) )
    stop("control$tol must be numeric")
  if ( !is.numeric(control$eps) )
    stop("control$eps must be numeric")
  if ( !is.logical(control$safe) )
    stop("control$eps must be logical")

  control
}

check.loaded.symbol <- function(symbol, dll="") {
  if ( !is.loaded(symbol, dll) ) {
    msg <- sprintf("Can't find C function %s", symbol)
    if ( dll != "" )
      sprintf("%s in shared library %s", msg, dll)
    stop(msg)
  }
  TRUE
}

check.info.ode <- function(info, control) {
  ## Might be doubling up here, but good to check.
  control <- check.control.ode(control)
  
  model <- if (is.null(info$name.ode)) info$name else info$name.ode
  if ( !(is.character(model) && length(model) == 1) )
    stop("'model' must be a single string")
  info$name.ode <- model
  info$ny    <- check.integer(check.scalar(info$ny))
  info$np    <- check.integer(check.scalar(info$np))
  info$idx.d <- check.integer(info$idx.d)

  if ( is.null(info$dll) )
    info$dll <- "diversitree"
  else if ( !(is.character(info$dll) && length(info$dll)) == 1 )
    stop("dll must be a single string")
  dll <- info$dll

  if ( control$compiled ) {
    if ( control$backend != "gslode" )
      stop("Only gslode backend supported with compiled code")
    check.loaded.symbol(sprintf("derivs_%s_gslode", model), dll)
  } else {
    if ( !is.function(info$derivs) )
      stop("info$derivs must be a function")
  }
  
  info
}

## For almost all models, there there must be a certain number of
## finite non-negative parameters.
check.pars.nonnegative <- function(pars, npar) {
  if ( length(pars) != npar )
    stop(sprintf("Incorrect parameter length: expected %d, got %d",
                 npar, length(pars)))
  if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
  pars
}

## Sometimes, it's useful just to check that the numbers look OK
## though:
check.nonnegative <- function(x, msg=NULL) {
  if ( is.null(msg) )
    msg <- "Parameters must be non-negative and finite"
  if ( any(!is.finite(x)) || any(x < 0)  )
    stop(msg)
  x
}
check.nonpositive <- function(x, msg=NULL) {
  if ( is.null(msg) )
    msg <- "Parameters must be non-positive and finite"
  if ( any(!is.finite(x)) || any(x > 0)  )
    stop(msg)
  x
}

## Check that a pointer is not NULL.
check.ptr <- function(ptr)
  .Call("check_ptr_not_null", ptr)
