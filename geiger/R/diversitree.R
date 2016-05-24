
## FUNCTIONS from diversitree 0.9-4: commit ec08b1d2bd from 06-22-2012 @ 'https://github.com/richfitz/diversitree'

.drop.tip=function(phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(phy)){
  
  
  if(missing(tip)) return(phy)
  if (is.character(tip)) tip <- which(phy$tip.label %in% tip)
  if(!length(tip)) return(phy)
    
  phy=as.phylo(phy)
  Ntip <- length(phy$tip.label)
  tip=tip[tip%in%c(1:Ntip)]
  if(!length(tip)) return(phy)


  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- nrow(phy$edge)

  wbl <- !is.null(phy$edge.length)
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !(edge2 %in% tip)  

  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel)) 
      break
    keep[sel] <- FALSE
  }

  phy2 <- phy
  phy2$edge <- phy2$edge[keep, ]
  if (wbl) 
    phy2$edge.length <- phy2$edge.length[keep]
  TERMS <- !(phy2$edge[, 2] %in% phy2$edge[, 1])
  oldNo.ofNewTips <- phy2$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  idx.old <- phy2$edge[TERMS, 2]
  phy2$edge[TERMS, 2] <- rank(phy2$edge[TERMS, 2])
  phy2$tip.label <- phy2$tip.label[-tip]
  if (!is.null(phy2$node.label))
    phy2$node.label <-
      phy2$node.label[sort(unique(phy2$edge[, 1])) - Ntip]
  phy2$Nnode <- nrow(phy2$edge) - n + 1L
  i <- phy2$edge > n
  phy2$edge[i] <- match(phy2$edge[i], sort(unique(phy2$edge[i]))) + n
  storage.mode(phy2$edge) <- "integer"
  collapse.singles(phy2)
}

argn <- function(x, ...)
UseMethod("argn")

`argn<-` <- function(x, value)
UseMethod("argn<-")

argn.bm <-
function(x, ...) {
	ret <- attr(x, "argn")
	if ( is.null(ret) )
    "beta" # Harmon
	else
    ret
}

.default.argn.mk2 <- function() c("q01", "q10")

.default.argn.mkn <- function(k) {
	base <- ceiling(log10(k + .5))
	fmt <- sprintf("q%%0%dd%%0%dd", base, base)
	sprintf(fmt, rep(1:k, each=k-1),
			unlist(lapply(1:k, function(i) (1:k)[-i])))
}

.format.levels.print=function(levels){
    base <- ceiling(log10(levels + .5))
	fmt <- sprintf("%%0%dd", base)
    fmt
}

.default.argn.bm <- function() "sigsq"


ROOT.FLAT  <- 1
ROOT.EQUI  <- 2
ROOT.OBS   <- 3
ROOT.GIVEN <- 4
ROOT.BOTH  <- 5
ROOT.MAX   <- 6
ROOT.ALL   <- ROOT.BOTH


.check.states <- function(phy, dat, strict.vals=NULL) {
	tree=phy
	states=dat
	allow.unnamed=FALSE
	strict=FALSE
	
	if ( is.matrix(states) ) {
## Multistate characters (experimental).  This will not work with
## clade trees, but they are only interesting for BiSSE, which has
## NA values for multistate (even weight).
		
		
		
		if ( inherits(tree, "clade.tree") )
		stop("Clade trees won't work with multistate tips yet")
		n <- rowSums(states > 0)
		if ( any(n == 0) )
		stop(sprintf("No state found for taxa: %s",
					 paste(names(tmp)[n == 0], collapse=", ")))
		
		i.mono <- which(n == 1)
		i.mult <- which(n >  1)
		
		tmp <- .matrix.to.list(states)
		names(tmp) <- rownames(states)
		
		states.mult <- lapply(tmp[i.mult], as.numeric)
		
		states <- rep(NA, length(tmp))
		names(states) <- names(tmp)
		states[i.mono] <- sapply(tmp[i.mono], function(x) which(x != 0))
		
		attr(states, "multistate") <- list(i=i.mult, states=states.mult)
	}
	
	if ( is.null(names(states)) ) {
		if ( allow.unnamed ) {
			if ( length(states) == length(tree$tip.label) ) {
				names(states) <- tree$tip.label
				warning("Assuming states are in tree$tip.label order")
			} else {
				stop(sprintf("Invalid states length (expected %s)", length(tree$tip.label)))
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
		if ( strict ) {
			if ( !isTRUE(all.equal(sort(strict.vals),
								   sort(unique(na.omit(states))))) )
			stop("Because strict state checking requested, all (and only) ",
				 sprintf("states in %s are allowed", paste(strict.vals, collapse=", ")))
		} else {
			extra <- setdiff(sort(unique(na.omit(states))), strict.vals)
			if ( length(extra) > 0 )
			stop(sprintf("Unknown states %s not allowed in states vector", paste(extra, collapse=", ")))
		}
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

.check.states.quasse <- function(tree, states, states.sd) {
	states <- .check.states(tree, states)
	
	if ( length(states.sd) == 1 )
		states.sd <- structure(rep(states.sd, length(states)), names=names(states))
	else
		states.sd <- .check.states(tree, states.sd)
	
	list(states=states, states.sd=states.sd)
}  



## Check that a number can reasonably be considered an integer.
.check.integer <- function(x) {
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


.check.tree <- function(tree, ultrametric=TRUE, bifurcating=TRUE, node.labels=FALSE) {
	if ( !inherits(tree, "phylo") )
    stop("'tree' must be a valid phylo tree")
	if ( ultrametric && !is.ultrametric(tree) )
    stop("'tree' must be ultrametric")
	if ( any(tree$edge.length < 0) )
    stop("Negative branch lengths in tree")
	## ape is.binary.tree() can let a few nasties through - for
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


.check.control.mkn <- function(control, k) {
	control <- modifyList(list(method="exp"), control)
	if ( control$method == "mk2" && k != 2 )
    stop("Method 'mk2' only valid when k=2")
	methods <- c("exp", "mk2", "ode")
	if ( !(control$method %in% methods) )
    stop(sprintf("control$method must be in %s",
                 paste(methods, collapse=", ")))
	control
}


.check.pars.nonnegative <- function(pars, npar) {
	if ( length(pars) != npar )
    stop(sprintf("Incorrect parameter length: expected %d, got %d",
                 npar, length(pars)))
	if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
	pars
}

.check.control.ode <- function(control=list()) {
	defaults <- list(safe=FALSE, tol=1e-8, eps=0, backend="deSolve", unsafe=FALSE)
	control <- modifyList(defaults, control)
	
#	backends <- c("deSolve", "cvodes", "CVODES")
	backends <- "deSolve"
	if ( length(control$backend) != 1 ) stop("'backend' must be a single element")
	control$backend <- backends[pmatch(control$backend, backends)]
	if ( is.na(control$backend) ) stop("Invalid backend selected")
	
#	if ( tolower(control$backend) == "cvodes" ) check.cvodes(error=TRUE)
	
	control$tol <- .check.scalar(control$tol)
	control$eps <- .check.scalar(control$eps)
	control$safe <- .check.scalar(control$safe)
	
	if ( !is.numeric(control$tol) )
    stop("control$tol must be numeric")
	if ( !is.numeric(control$eps) )
    stop("control$eps must be numeric")
	if ( !is.logical(control$safe) )
    stop("control$eps must be logical")
	
	control
}


.check.scalar <- function(x) {
	if ( length(x) != 1 )
    stop(deparse(substitute(x)), " must be a scalar")
	x
}

.check.info.ode <- function(info, control=NULL) {
	model <- if (is.null(info$name.ode)) info$name else info$name.ode
	if ( !(is.character(model) && length(model) == 1) )
    stop("'model' must be a single string")
	info$name.ode <- model
	info$ny    <- .check.integer(.check.scalar(info$ny))
	info$np    <- .check.integer(.check.scalar(info$np))
	info$idx.d <- .check.integer(info$idx.d)
	
	if ( is.null(info$dll) )
    info$dll <- "geiger"
	else if ( !(is.character(info$dll) && length(info$dll)) == 1 )
    stop("dll must be a single string")
	dll <- info$dll
	if ( is.null(info$banded) )
    info$banded <- FALSE
	else if ( length(info$banded) != 1 || !is.logical(info$banded) )
    stop("Invalid value for banded")
	
	if ( !is.null(control) ) {
		backend <- control$backend
		if ( is.function(info$derivs) ) {
			backend <- "R" # no effect.
		} else if ( backend == "deSolve" ) {
			.check.loaded.symbol(sprintf("initmod_%s", model), dll)
			.check.loaded.symbol(sprintf("derivs_%s", model),  dll)
		} else if ( backend == "cvodes" ) {
			.check.loaded.symbol(sprintf("derivs_%s_cvode", model), dll)
		} else if ( backend == "CVODES" ) {
			.check.loaded.symbol(sprintf("derivs_%s_cvode", model),       dll)
			.check.loaded.symbol(sprintf("initial_conditions_%s", model), dll)
		}
	}
	
	
	info
}

.check.loaded.symbol <- function(symbol, dll="") {
	if ( !is.loaded(symbol, dll) ) {
		msg <- sprintf("Cannot locate C function %s", symbol)
		if ( dll != "" )
		sprintf("%s in shared library %s", msg, dll)
		stop(msg)
	}
	TRUE
}


.make.cache <- function(tree) {
	## This works around some ape bugs with class inheritance.
	if (inherits(tree, "phylo")) class(tree) <- "phylo"
	edge <- tree$edge
	edge.length <- tree$edge.length
	idx <- seq_len(max(edge))
	n.tip <- length(tree$tip.label)
	tips <- seq_len(n.tip)
	root <- n.tip + 1
	
	is.tip <- idx <= n.tip
	
	children <- .get.children(edge, n.tip)
	
	parent <- edge[match(idx, edge[,2]),1]
	
	order <- .get.ordering(children, is.tip, root)
	len <- edge.length[match(idx, edge[,2])]
	
	## This is a bit of a hack, but this is to ensure that we can really
	## compute the depths accurately - this is a problem when there
	## joins (under split models) that occur right around nodes.
	height <- .branching.heights(tree)
	depth <- max(height) - height
	depth2 <- .branching.depth(len, children, order, tips)
	i <- abs(depth - depth2) < 1e-8
	depth[i] <- depth2[i]
	
	if ( is.ultrametric(tree) )
	## It is possible that an ultrametric tree will not quite have the
	## tips around zero.  This ensures it, which is is required for
	## dt.tips.grouped to work at present.
    depth[tips] <- 0
	
	## TODO: I do not need this ancestor thing for much - drop it here
	## and move it to the asr code that actually uses it (this takes a
	## lot of time, and is only used by the ASR code).
	## The only place that this is used at all is do.asr.marginal(); it
	## would be possible to make this as needed when making an
	## asr.marginal() function.
#	anc <- vector("list", max(order))
#	for ( i in c(rev(order[-length(order)]), tips) )
#   anc[[i]] <- c(parent[i], anc[[parent[i]]])
	
	ans <- list(tip.label=tree$tip.label,
				node.label=tree$node.label,
				len=len,
				children=children,
				parent=parent,
				order=order,
				root=root,
				n.tip=n.tip,
				n.node=tree$Nnode,
				tips=tips,
				height=height,
				depth=depth,
#				ancestors=anc,
				edge=edge,
				edge.length=edge.length)
	ans
}



.make.cache.bm <- function(tree, states, states.sd, control) {
	method <- control$method
	tree <- .check.tree(tree, ultrametric=FALSE)
	
	cache <- .make.cache(tree)
	if ( is.null(states.sd) || all(states.sd == 0) ) {
		cache$states <- .check.states(tree, states)
		cache$states.sd <- rep(0, cache$n.tip)
	} else {
		tmp <- .check.states.quasse(tree, states, states.sd)
		cache$states    <- tmp$states
		cache$states.sd <- tmp$states.sd
	}
	
	if ( method == "vcv" )
    cache$vcv <- vcv.phylo(tree)
	else
    cache$y <- .initial.tip.bm.pruning(cache)
	cache$info <- .make.info.bm(tree)
	cache
}

.make.info.bm <- function(phy) {
	list(name="bm",
		 name.pretty="Brownian motion",
## Parameters:
		 np=1L,
		 argn=.default.argn.bm(),
## Variables:
		 ny=3L,
		 k=NA,
		 idx.e=NA,
		 idx.d=NA,
## Phylogeny:
		 phy=phy,
## Inference:
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE,
## These are optional
		 doc=NULL,
		 reference=c(
			"Felsenstein, J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. American Journal of Human Genetics 25:471-492.")
	)
}



.make.cache.mkn <- function(tree, states, k, strict, control, ...) {
	method <- control$method
	method=match.arg(method, "exp")
	tree <- .check.tree(tree, ...)
	if ( !is.null(states) ) # for multitrait
    states <- .check.states(tree, states, strict.vals=1:k)
	cache <- .make.cache(tree)
#	if ( method == "mk2" )
#   cache$info <- .make.info.mk2(tree)
#	else
    cache$info <- .make.info.mkn(k, tree)
	cache$states  <- states
#	if ( method == "ode" ) {
#		cache$y <- initial.tip.mkn.ode(cache)
#		cache$info$name.ode <- "mknode"
#	}
	
	cache
}



.make.info.mkn <- function(k, phy) {
	list(name="mkn",
		 name.pretty="Mk(n)",
	## Parameters:
		 np=as.integer(k * k),
		 argn=.default.argn.mkn(k),
	## Variables:
		 ny=as.integer(k), # TODO/NEW: only for ode version...
		 k=as.integer(k),
		 idx.e=integer(0),
		 idx.d=seq_len(k),
	## Phylogeny:
		 phy=phy,
	## Inference:
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE,
	## These are optional
		 doc=NULL,
		 reference=c(
					 "Pagel (1994)",
					 "Lewis (2001)"))
}

.make.pars.mkn <- function(k) {
	qmat <- matrix(0, k, k)
	idx <- cbind(rep(1:k, each=k-1),
				 unlist(lapply(1:k, function(i) (1:k)[-i])))
	npar <- k*(k-1)
	
	function(pars) {
		.check.pars.nonnegative(pars, npar)
		qmat[idx] <- pars
		diag(qmat) <- -rowSums(qmat)
		qmat
	}
}


.make.pij.mkn <- function(info, control) {
	method=match.arg(control$method, "exp")
	control <- .check.control.ode(control)
	## TODO/NEW: cvodes should be easy too.
	if ( control$backend != "deSolve" )
    stop("Only deSolve backend is available")
	
	k <- info$k
	info <- list(name="mkn_pij",ny=k*k, np=k*k, idx.d=integer(0))
	info <- .check.info.ode(info, control)
	pij.ode <- .make.ode.deSolve(info, control)
	
	yi <- diag(k) # initial conditions always same.
	
	function(len, pars)
    pij.ode(yi, c(0, len), pars)
}

.make.ode.deSolve <- function(info, control) {
	model  <- info$name.ode
	np     <- info$np  
	ny     <- info$ny
	dll    <- info$dll
	banded <- info$banded
	
	safe   <- control$safe
	unsafe <- control$unsafe
	atol   <- rtol <- as.numeric(control$tol)
	
	initfunc <- sprintf("initmod_%s", model)
	derivs <- sprintf("derivs_%s", model)
	
    safe=TRUE # to eliminate foreign-function call for CRAN check ## JME 03.27.2013
    
	if ( safe ) {
		function(vars, times, pars) {
			ret <- t(lsoda(vars, times, pars, rtol=rtol, atol=atol,
						   initfunc=initfunc, func=derivs,
						   dllname=dll)[-1,-1,drop=FALSE])
			dimnames(ret) <- NULL
			ret
		}
	} else {
## Temporary fix so that I can work on the cluster.  This will be
## removed and DESCRIPTION updated to require R 2.12.0 or greater.
		if ( getRversion() >= package_version("2.12.0") )
		vers <-  packageVersion("deSolve")
		else
		vers <- package_version(packageDescription("deSolve",
												   fields="Version"))
		max.deSolve <- package_version("1.10-4")
		if ( !unsafe && vers > max.deSolve ) {
			str <- paste(sQuote("geiger"), "is not known to work with deSolve >",
						 max.deSolve, "\n\tfalling back on slow version")
			warning(str)
			control$safe <- TRUE
			return(.make.ode.deSolve(info, control))
		}
		
		initfunc  <- getNativeSymbolInfo(initfunc, PACKAGE=dll)$address
		derivs    <- getNativeSymbolInfo(derivs,   PACKAGE=dll)$address
		jacfunc <- NULL
## Magic numbers from deSolve:lsoda.R
		jactype <- if ( banded ) 5L else 2L
		
		maxordn <- 12
		maxords <- 5
		lrw <- as.integer(max(20 + ny * (maxordn + 1) + 3 * ny,
							  20 + ny * (maxords + 1) + 3 * ny + ny^2 + 2))
		liw <- as.integer(20 + ny)
		
		iwork <- vector("integer", 20)
		iwork[1] <- as.integer(1)    # banddown
		iwork[2] <- as.integer(1)    # bandup
		iwork[6] <- as.integer(5000) # maxsteps
		
		rwork <- vector("double", 20)
		rwork[5] <- 0                     # hini
		rwork[6] <- 10                    # hmax (consider 0)
		rwork[7] <- 0                     # hmin
		
		INTZERO <- 0L
		INTONE <- 1L
		INTTWO <- 2L
		
		flist <- list(fmat=0, tmat=0, imat=0, ModelForc=NULL)
		elag <- list(islag=0L)
		
##----- eliminated foreign-function call for CRAN check ## JME 03.27.2013
#		sol <- function(vars, times, pars) {
#			if ( length(vars) != ny )
#			stop("Incorrect variable length")
#			if ( length(times) <= 1 )
#			stop("Need >= 2 times")
#			storage.mode(vars) <- storage.mode(times) <- "numeric"

#			ret <- 
#			.Call("call_lsoda", vars, times, derivs, pars,
#				  rtol, atol,
#				  NULL,      # rho: environment
#				  NULL,      # tcrit: critical times
#				  jacfunc,
#				  initfunc,
#				  NULL,      # eventfunc [New in 1.6]
#				  INTZERO,   # verbose (false)
#				  INTONE,    # itask
#				  rwork,
#				  iwork,
#				  jactype,   # Jacobian type (2=full, 5=banded [1 up and down])
#				  INTZERO,   # nOut (no global variables)
#				  lrw,       # size of workspace (real)
#				  liw,       # size of workspace (int)
#				  INTONE,    # Solver
#				  NULL,      # rootfunc
#				  INTZERO,   # nRoot
#				  0,         # rpar: no extra real parameters
#				  INTZERO,   # ipar: no extra integer parameters
#				  INTZERO,   # Type
#				  flist,     # [New in 1.5]
#				  list(),    # elist [New in 1.6]
#				  elag,      # [New in 1.7]
#				  PACKAGE="deSolve")
#			if ( max(abs(ret[1,] - times)) > 1e-6 )
#			stop("Integration error: integration stopped prematurely")
#			ret[-1,-1,drop=FALSE]
#		}
    stop("unsupported call")
	}
}


.root.bm.direct <- function(vars, log.comp, root, root.x) {
	if ( root == "max" ) {
## This treats a prior on the root as a delta function centred at
## the ML root state.
## The first term can be more intuitively written as:
##   dnorm(vars[1], vars[1], sqrt(vars[2]), TRUE)
##   dnorm(0, 0, sqrt(vars[2]), TRUE)
		- log(2 * pi * vars[[2]]) / 2 + vars[[3]] + sum(log.comp)
	} else if ( root == "flat" ) {
## Flat prior (by this point, function integrates to vars[[3]])
		vars[[3]] + sum(log.comp)
	} else if ( root == "obs" ) {
## Observed weighting (integrate norm norm wrt x from -inf to inf
## gives 1 / (2 sqrt(pi s2))).
		-log(2 * sqrt(pi * vars[[2]])) + vars[[3]] + sum(log.comp)
	} else if ( root == "given" ) {
		dnorm(root.x, vars[1], sqrt(vars[2]), TRUE) + vars[[3]] +
		sum(log.comp)
	} else {
		stop("Invalid root mode")
	}
}


.dt.tips.ordered <- function(y, tips, t) {
	i <- order(t)
	
	if ( is.list(y) ) {
		if ( length(y) != length(tips) )
		stop("y must be same length as tips")
		if ( length(y) != length(t) )
		stop("y must be the same length as t")
		list(target=tips[i],
			 t=t[i],
			 y=y[i])
	} else if ( is.matrix(y) ) {
		if ( ncol(y) != length(tips) )
		stop("y must be same length as tips")
		if ( ncol(y) != length(t) )
		stop("y must be the same length as t")
		list(target=tips[i],
			 t=t[i],
			 y=y[,i], type="ORDERED")
	} else {
		stop("y must be a list or matrix")
	}
}




is.constrained <- function(x) inherits(x, "constrained")

## First up, consider the one-shot case: do not worry about incremental
## updates.
## 
## For the first case, everything is OK on the lhs and rhs
## For subsequent cases:
##   lhs cannot contain things that are
##      - constrained things (already lhs anywhere)
##      - things constrained to (things on the rhs anywhere)
##   rhs cannot contain things that are
##      - constrained things (already lhs anywhere)
## It is possibly worth pulling out all the numerical constants and
## the "paired" parameters here to avoid using eval where
## unnecessary.  However, this makes the function substantially uglier
## for a very minor speedup.
.constrain <- function(f, ..., formulae=NULL, names=argn(f), extra=NULL) {
	if ( is.constrained(f) ) {
        stop("'f' appears already constrained")
		formulae <- c(attr(f, "formulae"), formulae)
		f <- attr(f, "func")
	}
    
    levels=attr(f, "levels")
	
	formulae <- c(formulae, list(...))
	names.lhs <- names.rhs <- names
	rels <- list()
	
	for ( formula in formulae ) {
		res <- .constrain.parse(formula, names.lhs, names.rhs, extra)
		if ( attr(res, "lhs.is.target") ) {
			i <- which(sapply(rels, function(x) identical(x, res[[1]])))
			rels[i] <- res[[2]]
			
## This will not work with *expressions* involving the LHS; that
## would require rewriting the expressions themselves (which
## would not be too hard to do).  But for now let us just cause
## an error...
			lhs.txt <- as.character(res[[1]])
			if ( any(sapply(rels, function(x) lhs.txt %in% all.vars(x))) )
			stop(sprintf("lhs (%s) is in an expression and can't be constrained",
						 lhs.txt))
		}
		
		names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
		names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
		rels <- c(rels, structure(res[2], names=as.character(res[[1]])))
	}
	
	i <- match(unique(sapply(rels, as.character)), extra)
	final <- c(extra[sort(i[!is.na(i)])], names.rhs)
	npar <- length(final)
	
## "free" are the parameters that have nothing special on their RHS
## and are therefore passed directly through
	free <- setdiff(names.rhs, names(rels))
	free.i <- match(free, names) # index in full variables
	free.j <- match(free, final) # index in given variables.
	
## Targets are processed in the same order as given by formulae. 
	target.i <- match(names(rels), names)
	
	pars.out <- rep(NA, length(names))
	names(pars.out) <- names
	g <- function(pars, ..., pars.only=FALSE) {
		if ( length(pars) != npar )
		stop(sprintf("Incorrect parameter length: expected %d, got %d",
					 npar, length(pars)))
		
		pars.out[free.i] <- pars[free.j]
		e <- structure(as.list(pars), names=final)
		pars.out[target.i] <- unlist(lapply(rels, eval, e))
		
		if ( pars.only )
		pars.out
		else
		f(pars.out, ...)
	}
	
	class(g) <- c("constrained", class(f))
	attr(g, "argn") <- final
	attr(g, "formulae") <- formulae
	attr(g, "extra") <- extra
	attr(g, "func") <- f
    attr(g, "levels") <- levels ## JME
	g
}


## Parsing constraints:
## The LHS of a formula must be a single variable name that exists in
## "names.lhs"
##
## The RHS can be one of
##   - numeric value
##   - expression
## If it is an expression, then all variable names must be found in
## names.rhs (or perhaps in the containing environment - check in the
## future?)

## I might eventually allow formulae of the form
##   lambda | lambda1 ~ lambda0
## to allow renaming?

## How do I use this to allow recasting to alternative bases?
## Forward definitions for just the diversification rate are
##   foo(f, r0 ~ lambda0 - mu0, r1 ~ lambda1 - mu1)
## and for both diversification and relative extinction    
##   foo(f,
##       r0 ~ lambda0 - mu0, r1 ~ lambda1 - mu1,
##       e0 ~ mu0 / lambda0, e1 ~ mu1 / lambda1)
## Backward definitions (leaving mu unchanged)
##   foo(f, lambda0 ~ r0 + mu0, lambda1 ~ r1 + mu1)
## both:
##   foo(f, lambda0 ~ r0/(1 - e0), lambda1 ~ r1/(1 - e1),
##       r0 * e0 / (1 - e0), r1 * e1 / (1 - e1))
.constrain.parse <- function(formula, names.lhs, names.rhs, extra=NULL) {
	formula <- as.formula(formula)
	if ( length(formula) != 3L )
    stop("Invalid formula")
	lhs <- formula[[2]]
	rhs <- formula[[3]]
	
## Checking the lhs is easy: is the lhs in the list of allowable
## names and of length 1?  Anything that does not match this is
## invalid.
	if ( !is.name(lhs) )
    stop("Invalid target on LHS of formula" )
	lhs.is.target <- is.na(match(as.character(lhs), names.lhs))
	
## Checking the rhs is more difficult.  We are OK if any of the
## following is met:
##   Numeric values (checked at the end)
##   If all vars in rhs are in names.rhs
##   There is a single var and it is in extra
## Failing that, if the rhs is a single variable that does exist in
## the calling environment.
	if ( is.language(rhs) ) {
		vars <- all.vars(rhs)
		ok <- (all(vars %in% names.rhs) ||
			   length(vars) == 1 && vars %in% extra)
		if ( !ok && length(vars) == 1 ) {
			e <- parent.frame()
			if ( exists(vars, e) ) {
				rhs <- get(vars, e)
				ok <- TRUE
			}
		}
		
		if ( !ok )
		stop("Invalid RHS of formula:\n\t", as.character(rhs))
		if ( as.character(lhs) %in% vars )
		stop("LHS cannot appear in RHS")
	} else if ( !is.numeric(rhs) ) {
		stop("RHS must be expression, variable or number")
	}
	res <- list(lhs, rhs)
	attr(res, "lhs.is.target") <- lhs.is.target
	res
}

.initial.tip.bm.pruning <- function(cache) {
	y <- mapply(function(mean, sd) c(mean, sd*sd, 0),
				cache$states, cache$states.sd, SIMPLIFY=TRUE)
	.dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
}



.rootfunc.mkn <- function(res, pars, root, root.p, intermediates) {
	d.root <- res$vals
	lq <- res$lq
	k <- length(d.root)
	
	root.p <- .root.p.calc(d.root, pars, root, root.p, NULL)
	if ( root == "all" )
    loglik <- log(d.root) + sum(lq)
	else
    loglik <- log(sum(root.p * d.root)) + sum(lq)
	
	if ( intermediates ) {
		res$root.p <- root.p
		attr(loglik, "intermediates") <- res
		attr(loglik, "vals") <- d.root
	}
	
	loglik
}


.root.p.calc <- function(vals, pars, root, root.p=NULL, root.equi=NULL) {
	if ( !is.null(root.p) && root != "given" )
    warning("Ignoring specified root state")
	
	k <- length(vals)
	
	if ( root == "flat" ) {
		p <- 1/k
	} else if ( root == "equi" ) {
		if ( is.null(root.equi) )
		stop("Equilibrium root probability not possible with this method")
		p <- root.equi(pars)
	} else if ( root == "obs") {
		p <- vals / sum(vals)
	} else if ( root == "given" ) {
		if ( length(root.p) != length(vals) )
		stop("Invalid length for root.p")
		p <- root.p
	} else if ( root == "all" ) {
		p <- rep(1, k)
	} else {
		stop("Invalid root mode")
	}
	p
}

.branching.heights <- function(phy) {
	if (!inherits(phy, "phylo"))
    stop('object "phy" is not of class "phylo"')
	
	edge <- phy$edge
	n.node <- phy$Nnode
	n.tip <- length(phy$tip.label)
	
	ht <- numeric(n.node + n.tip) # zero
	for (i in seq_len(nrow(edge)) )
    ht[edge[i, 2]] <- ht[edge[i, 1]] + phy$edge.length[i]
	
## Ugly, but fairly compatible with branching.times()
	names.node <- phy$node.label
	if (is.null(names.node))
    names.node <- (n.tip + 1):(n.tip + n.node)
	names(ht) <- c(phy$tip.label, names.node)
	
	ht
}

.branching.depth <- function(len, children, order, tips) {
	depth <- numeric(nrow(children))
	depth[tips] <- 0
	for ( i in order )
    depth[i] <- depth[children[i,1]] + len[children[i,1]]
	depth
}

.get.children <- function(edge, n.tip) {
## To construct the children vector, this is what I used to do:
##   lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
## But this is slow for large trees.  This is faster:
## children <- split(edge[,2], edge[,1])
## Surprisingly, most of the time is in coercing edge[,1] into a
## factor.
	x <- as.integer(edge[,1])
	levels <- as.integer((n.tip+1):max(edge[,1]))
	f <- match(x, levels)
	levels(f) <- as.character(levels)
	class(f) <- "factor"
	children <- split(edge[,2], f)
	names(children) <- NULL
	
## In most cases, this will have been checked by .check.tree()
## This is currently the time sink here.
	if ( !all(unlist(lapply(children, length)) == 2) )
    stop("Multifircations/unbranched nodes in tree - best get rid of them")
	
	rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), 2)))
}

.get.ordering <- function(children, is.tip, root) {
	todo <- list(root)
	i <- root
	repeat {
		kids <- children[i,]
		i <- kids[!is.tip[kids]]
		if ( length(i) > 0 )
		todo <- c(todo, list(i))
		else
		break
	}
	as.vector(unlist(rev(todo)))
}


.matrix.to.list <- function(m) {
	n <- nrow(m)
	out <- vector("list", n)
	for ( i in seq_len(n) )
    out[[i]] <- m[i,]
	out
}


.toC.int <- function(x) {
	x <- x - 1
	storage.mode(x) <- "integer"
	x
}



