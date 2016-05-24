# Development of a dorng equivalent to dopar for reproducible loops
# 
# Author: Renaud Gaujoux
# Creation: 17 Aug 2011
###############################################################################

#library(foreach)

# or-NULL operator (borrowed from Hadley Wickham)
'%||%' <- function(x, y) if( !is.null(x) ) x else y

.collapse <- function(x, n=length(x), sep=', '){
	
	res <- paste(if( missing(n) ) x else head(x, n), collapse=', ')
	if( length(x) > n )
		res <- paste(res, '...', sep=', ')
	res
}

#' Back Compatibility Option for doRNG
#' 
#' Sets the behaviour of \%dorng\% foreach loops from a
#' given version number.
#' 
#' @section Behaviour changes in versions:
#' 
#' \describe{
#' \item{1.4}{ The behaviour of \code{doRNGseed}, and therefore of 
#' \code{\%dorng\%} loops, changed in the case where the current RNG was 
#' L'Ecuyer-CMRG.
#' Using \code{set.seed} before a non-seeded loop used not to be identical 
#' to seeding via \code{.options.RNG}.
#' Another bug was that non-seeded loops would share most of their RNG seed!
#' }
#' }
#' 
#' @param x version number to switch to, or missing to get the currently 
#' active version number, or \code{NULL} to reset to the default behaviour, 
#' i.e. of the latest version.
#' 
#' @return a character string
#' If \code{x} is missing this function returns the version number from the 
#' current behaviour.
#' If \code{x} is specified, the function returns the old value of the 
#' version number (invisible).
#' 
#' @importFrom utils packageVersion
#' @export
#' @examples 
#'
#' \dontshow{ registerDoSEQ() }
#' 
#' ## Seeding when current RNG is L'Ecuyer-CMRG
#' RNGkind("L'Ecuyer")
#' 
#' doRNGversion("1.4")
#' # in version >= 1.4 seeding behaviour changed to fix a bug
#' set.seed(123)
#' res <- foreach(i=1:3) %dorng% runif(1)
#' res2 <- foreach(i=1:3) %dorng% runif(1)
#' stopifnot( !identical(attr(res, 'rng')[2:3], attr(res2, 'rng')[1:2]) )
#' res3 <- foreach(i=1:3, .options.RNG=123) %dorng% runif(1)
#' stopifnot( identical(res, res3) )
#' 
#' # buggy behaviour in version < 1.4
#' doRNGversion("1.3")
#' res <- foreach(i=1:3) %dorng% runif(1)
#' res2 <- foreach(i=1:3) %dorng% runif(1)
#' stopifnot( identical(attr(res, 'rng')[2:3], attr(res2, 'rng')[1:2]) )
#' res3 <- foreach(i=1:3, .options.RNG=123) %dorng% runif(1)
#' stopifnot( !identical(res, res3) )
#' 
#' # restore default RNG  
#' RNGkind("default")
#' # restore to current doRNG version
#' doRNGversion(NULL)
#' 
doRNGversion <- local({

	currentV <- "1.5.3" #as.character(packageVersion('doRNG')) 
	cache <- currentV
	function(x){
		if( missing(x) ) return(cache)
		if( is.null(x) ) x <- currentV 
		
		# update cache and return old value
		old <- cache
		cache <<- x
		invisible(old)
	}
})

#' @importFrom utils compareVersion
checkRNGversion <- function(x){
	compareVersion(doRNGversion(), x)
}

doRNGseq <- function(n, seed=NULL, ...){
	
	# compute sequence using rngtools::RNGseq
#	library(rngtools)
	res <- RNGseq(n, seed, ..., version=if( checkRNGversion('1.4') >=0 ) 2 else 1)
	
}

#' Getting Information About doRNG Foreach Backend
#' 
#' \code{infoDoRNG} returns information about the doRNG backend, e.g., version,
#' number of workers.
#' It is not meant to be called by the user.
#' 
#' 
#' @param data a list of data used by the backend
#' @param item the data item requested, as a character string 
#' (e.g. 'name', 'workers', 'version')
#' 
#' @return \code{infoDoRNG} returns the requested info (usually as a character 
#' string or a numeric value).
#' 
#' @keywords internal
#' @author Renaud Gaujoux
#' 
infoDoRNG <- function (data, item) 
{	
	switch(item
			, workers = data$backend$info(data$backend$data, "workers")
			, name = "doRNG"
			, version = "doRNG 1.5.3" 
			, NULL)
}

#' \code{doRNG} implements the generic reproducible foreach backend. It should 
#' not be called directly by the user.
#' 
#' @param obj a foreach description of the loop arguments
#' @param ex the lopp expression
#' @param envir the loop's evaluation environment
#' @param data configuration data of the doRNG backend
#' 
#' @rdname infoDoRNG
doRNG <- function (obj, ex, envir, data){
		
			
	if( is.null(obj$options) )
		obj$options <- list()
	
	if( !'RNG' %in% names(obj$options) ){
		obj$options$RNG <- if( !data$once || data$nseed==0 ){
			#message("doRNG backend - use seed ", if( data$once ) "only once" else "for every loop", ":")
			data$seed
		}
		else
			NULL
	}
	
#	data$nseed <- data$nseed + 1					
#	assign('data', data, pos=foreach:::.foreachGlobals)	
	
	rngBackend <- getDoBackend()
    # increment number of calls to doRNG
    rngBackend$data$nseed <- rngBackend$data$nseed + 1
    # directly register (temporarly) the computing backend
	on.exit({setDoBackend(rngBackend)}, add=TRUE)
	setDoBackend(rngBackend$data$backend)
	do.call('%dorng%', list(obj, ex), envir=parent.frame())
}

##% Get/Sets the registered foreach backend's data
getDoBackend <- function(){
    # one has to get the complete set of backend data from within the foreach Namespace
    foreach_ns <- asNamespace('foreach')
#    .foreachGlobals <- get('.foreachGlobals', foreach_ns)
    .foreachGlobals <- ns_get('.foreachGlobals', foreach_ns)
#    getDoPar <- get('getDoPar', foreach_ns)
    getDoPar <- ns_get('getDoPar', foreach_ns)
	c(getDoPar()
	, info= if( exists("info", where = .foreachGlobals, inherits = FALSE) ) .foreachGlobals$info else function(data, item) NULL)
}
setDoBackend <- function(backend){
	ob <- getDoBackend()
	do.call('setDoPar', backend)
	invisible(ob)
}

#' Reproducible Parallel Foreach Backend
#' 
#' \code{\%dorng\%} is a foreach operator that provides an alternative operator 
#' \code{\%dopar\%}, which enable reproducible foreach loops to be performed.
#' 
#' @param obj a foreach object as returned by a call to \code{\link{foreach}}.
#' @param ex the \code{R} expression to evaluate.
#' 
#' @return \code{\%dorng\%} returns the result of the foreach loop. See \code{\link{\%dopar\%}}.
#' The whole sequence of RNG seeds is stored in the result object as an attribute.
#' Use \code{attr(res, 'rng')} to retrieve it. 
#' 
#' @importFrom iterators iter
#' @export 
#' @rdname doRNG
#' @usage obj \%dorng\% ex
#' @seealso \code{\link{foreach}}, \code{\link[doParallel]{doParallel}}
#' , \code{\link[doParallel]{registerDoParallel}}, \code{\link[doMPI]{doMPI}}
#' @examples 
#' 
#' library(doParallel)
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#' 
#' # standard %dopar% loops are _not_ reproducible
#' set.seed(1234)
#' s1 <- foreach(i=1:4) %dopar% { runif(1) }
#' set.seed(1234)
#' s2 <- foreach(i=1:4) %dopar% { runif(1) }
#' identical(s1, s2)
#' 
#' # single %dorng% loops are reproducible
#' r1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' r2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' identical(r1, r2)
#' # the sequence os RNG seed is stored as an attribute
#' attr(r1, 'rng')
#' 
#' # stop cluster
#' stopCluster(cl)
#' 
#' # More examples can be found in demo `doRNG`
#' \dontrun{
#' demo('doRNG')
#' }
#' 
#' @demo Some features of the %dorng% foreach operator 
#' 
#' library(doRNG)
#' library(doParallel)
#' registerDoParallel(cores=2)
#' 
#' # single %dorng% loops are reproducible
#' r1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' r2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' identical(r1, r2)
#' # the sequence os RNG seed is stored as an attribute
#' attr(r1, 'rng')
#' 
#' # sequences of %dorng% loops are reproducible
#' set.seed(1234)
#' s1 <- foreach(i=1:4) %dorng% { runif(1) }
#' s2 <- foreach(i=1:4) %dorng% { runif(1) }
#' # two consecutive (unseed) %dorng% loops are not identical
#' identical(s1, s2)
#'
#' # But the whole sequence of loops is reproducible
#' set.seed(1234)
#' s1.2 <- foreach(i=1:4) %dorng% { runif(1) }
#' s2.2 <- foreach(i=1:4) %dorng% { runif(1) }
#' identical(s1, s1.2) && identical(s2, s2.2)
#' # it gives the same result as with .options.RNG
#' identical(r1, s1) 
#' 
#' # Works with SNOW-like and MPI clusters
#' # SNOW-like cluster
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#' 
#' s1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' s2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' identical(s1, s2)
#' 
#' stopCluster(cl)
#' registerDoSEQ()
#' 
#' # MPI cluster
#' library(doMPI)
#' cl <- startMPIcluster(2)
#' registerDoMPI(cl)
#' 
#' s1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' s2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' identical(s1, s2)
#' 
#' closeCluster(cl)
#' registerDoSEQ()
#' 
#' 
`%dorng%` <- function(obj, ex){
	
	#library(rngtools)
	
	# exit if nested or conditional loop
	if( any(c('xforeach', 'filteredforeach') %in% class(obj)) )
		stop("nested/conditional foreach loops are not supported yet.\nSee the package's vignette for a work around.")
	
	# if an RNG seed is provided then setup random streams 
	# and add the list of RNGs to use as an iterated arguments for %dopar%
#	library(parallel)
	obj$argnames <- c(obj$argnames, '.doRNG.stream')
	it <- iter(obj)
	argList <- as.list(it)
	
	# restore current RNG  on exit if a seed is passed
	rngSeed <- 
	if( !is.null(obj$options$RNG) ){
		
		# setup current RNG restoration
		RNG.old <- RNGseed()
		on.exit({RNGseed(RNG.old)}, add=TRUE)
		
		# extract RNG setting from object if possible (do not resolve single seed)
		rngSeed <- getRNG(obj$options$RNG, num.ok=TRUE) %||% obj$options$RNG
		
		# ensure it is a list
		# NB: unnamed lists are sequences of seeds
		if( !is.list(rngSeed) || is.null(names(rngSeed)) ){
			rngSeed <- list(rngSeed)
		}
		rngSeed
	}
	
	# generate a sequence of streams
#	print("before RNGseq")
#	print(head(RNGseed()))
	obj$args$.doRNG.stream <- do.call("doRNGseq", c(list(n=length(argList), verbose=obj$verbose), rngSeed))
#	print("after RNGseq")
#	print(head(RNGseed()))
	#print(obj$args$.doRNG.stream)
	## SEPCIAL CASE FOR doSEQ or doMPI
	# TODO: figure out why doMPI draw once from the current RNG (must be linked
	# to using own code to setup L'Ecuyer RNG)
	# restore RNG settings as after RNGseq if doSEQ is the backend and no seed was passed
	dp <- getDoParName()
#	print(dp)
	if( is.null(obj$options$RNG) && (is.null(dp) || dp=='doSEQ' || dp=='doMPI') ){
#		print("reset as after RNGseq")
		RNG.old <- RNGseed()
		on.exit({RNGseed(RNG.old)}, add=TRUE)
	}
	##
	
	# export package doRNG if not already exported
	if( !('doRNG' %in% obj$packages) )
		obj$packages <- c(obj$packages, 'doRNG')
	
	# append code to the loop expression to set the RNG
	ex <- as.call(list(as.name('{'),
					quote({rngtools::RNGseed(.doRNG.stream);}),
					substitute(ex)))
	
	# directly register (temporarly) the computing backend
	if( !is.null(dp) && dp == 'doRNG' ){
		rngBackend <- getDoBackend()
		on.exit({setDoBackend(rngBackend)}, add=TRUE)
		setDoBackend(rngBackend$data$backend)
	}
	
	# call the standard %dopar% operator
	res <- do.call('%dopar%', list(obj, ex), envir=parent.frame())
	# add seed sequence as an attribute
	attr(res, 'rng') <- obj$args$.doRNG.stream
	# return result
	res
}

#' Registering doRNG for Persistent Reproducible Parallel Foreach Loops   
#' 
#' \code{registerDoRNG} registers the doRNG foreach backend.
#' Subsequent \%dopar\% loops are then performed using the previously 
#' registered foreach backend, but are internally performed as \code{\link{\%dorng\%}} loops,
#' making them fully reproducible.
#' 
#' Briefly, the RNG is set, before each iteration, with seeds for L'Ecuyer's CMRG 
#' that overall generate a reproducible sequence of statistically independent 
#' random streams.
#' 
#' Note that (re-)registering a foreach backend other than doRNG, after a call 
#' to \code{registerDoRNG} disables doRNG -- which then needs to be registered.
#' 
#' @param seed a numerical seed to use (as a single or 6-length numerical value)
#' @param once a logical to indicate if the RNG sequence should be seeded at the 
#' beginning of each loop or only at the first loop. 
#' 
#' @seealso \code{\link{\%dorng\%}}
#' @export
#' @examples 
#' 
#' library(doParallel)
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#' 
#' # One can make reproducible loops using the %dorng% operator
#' r1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
#' # or convert %dopar% loops using registerDoRNG
#' registerDoRNG(1234)
#' r2 <- foreach(i=1:4) %dopar% { runif(1) }
#' identical(r1, r2)
#' stopCluster(cl)
#'
#' # Registering another foreach backend disables doRNG
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#' set.seed(1234)
#' s1 <- foreach(i=1:4) %dopar% { runif(1) }
#' set.seed(1234)
#' s2 <- foreach(i=1:4) %dopar% { runif(1) }
#' identical(s1, s2)
#' \dontshow{ stopifnot(!identical(s1, s2)) }
#' 
#' # doRNG is re-nabled by re-registering it 
#' registerDoRNG()
#' set.seed(1234)
#' r3 <- foreach(i=1:4) %dopar% { runif(1) }
#' identical(r2, r3)
#' # NB: the results are identical independently of the task scheduling
#' # (r2 used 2 nodes, while r3 used 3 nodes)
#' 
#' # argument `once=FALSE` reseeds doRNG's seed at the beginning of each loop 
#' registerDoRNG(1234, once=FALSE)
#' r1 <- foreach(i=1:4) %dopar% { runif(1) }
#' r2 <- foreach(i=1:4) %dopar% { runif(1) }
#' identical(r1, r2)
#' 
#' # Once doRNG is registered the seed can also be passed as an option to %dopar%
#' r1.2 <- foreach(i=1:4, .options.RNG=456) %dopar% { runif(1) }
#' r2.2 <- foreach(i=1:4, .options.RNG=456) %dopar% { runif(1) }
#' identical(r1.2, r2.2) && !identical(r1.2, r1)
#' \dontshow{ stopifnot(identical(r1.2, r2.2) && !identical(r1.2, r1)) }
#' 
#' stopCluster(cl)
#' 
registerDoRNG <- function(seed=NULL, once=TRUE){
	
	backend <- getDoBackend()
	
	# use stored backend if registerDoRNG was called repeatedly
	if( getDoParName() == 'doRNG' )
		backend <- backend$data$backend
	
	# set the current RNG with seed immediately if only used once
	if( once && !is.null(seed) ){
		if( !is.numeric(seed) || length(seed)!=1L )
			stop("Invalid seed: must be a single numeric value.")
		set.seed(seed)
		seed <- NULL
	}
	
	setDoPar(doRNG, list(seed=seed, once=once, nseed=0, backend=backend), infoDoRNG)	
	
}

###% Reproducibly Apply a Function over a List or Vector
###% 
###% @aliases xapply reproduce
###% 
###% \code{reproduce} and \code{xapply} are a reproducible versions 
###% of \code{\link{replicate}} and \code{\link{sapply}} respectively, 
###% that ensures the reproducibility of the results, when stochastic computations
###% are involved.
###% 
###% The reproducibility is achieved by using LEcuyer's RNG provided by R core
###% since R-2.14.0, to generate independent random streams 
###% that are used as the random number generator for each replicate.
###% 
###% @param n the number of replication as a single numeric (integer)
###% @param seed the main numerical seed used to initialize the sequence of random 
###% streams
###% @param expr the expression (language object, usually a call) to evaluate repeatedly
###% @param simplify logical; should the result be simplified to a vector or 
###% matrix if possible?
###% 
###% 
###%  
#reproduce <- function (n, expr, seed=NULL, simplify = TRUE){
#	f <- eval.parent(substitute(function(...) expr))
#	xapply(integer(n), seed, f, simplify = simplify)
#}
#
#xapply <- function (X, FUN, seed=NULL, ..., simplify = TRUE, USE.NAMES = TRUE){
#	
#	# generate a sequence of streams
#	.RNG.stream <- RNGseq(length(X), seed, packed=TRUE)
#	
#	# keep current RNG and restore it on exit (useful for the sequential backend doSEQ)
#	RNG.old <- rstream.RNG()
#	on.exit(rstream.RNG(RNG.old), add=TRUE)
#	
#	# append code to the loop expression to set the RNG	
#	expr <- as.call(list(as.name('{'),
#					quote({doRNGseed(.rng);}),
#					quote(do.call(FUN, list(...)))))
#	
#	env <- environment(FUN)
#	f <- eval(substitute(function(.rng, ..., FUN) expr), env)	
#	mapply(f, .RNG.stream, X, MoreArgs=c(list(...), FUN=FUN), 
#			SIMPLIFY = simplify, USE.NAMES= USE.NAMES)
#}
