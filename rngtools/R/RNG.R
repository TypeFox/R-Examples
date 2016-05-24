# Copyright (C) 2009-2012 Renaud Gaujoux
# 
# This file is part of the rngtools package for R. 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
# Creation: 08 Nov 2011
###############################################################################

library(pkgmaker)

###% Returns all the libraries that provides a user-supplied RNG
###% 
###% The library that provides the wrapper hooks for the management multiple 
###% user-supplied RNG is removed from the output list.
###% 
RNGlibs <- function(n=0, full=FALSE, hook="user_unif_rand", unlist=TRUE){
	dlls <- getLoadedDLLs()
	res <- lapply(dlls, function(d){
				dname <- d[['name']]
				if( dname=='' )
					return(NA)
				
				symb.unif_rand <- RNGlib(PACKAGE=dname, hook=hook)
				if( is.null(symb.unif_rand) )
					NA
				else
					symb.unif_rand
			})
	
	res <- res[!is.na(res)]
	if( !full )
		res <- names(res)	
	
	# limit the results if requested
	if( n>0 )
		res <- tail(res, n)
	
	# return result
	if( unlist && length(res) == 1 )
		res[[1]]
	else
		res
}

###% Returns the library that provides the current user-supplied RNG hooks.
###% 
###% This is the library that is first called by runif when using setting RNG 
###% kind to "user-supplied".
###% In general this will be rstream, except if a package providing the RNG hook 
###% 'user_unif_rand' is loaded after rstream, and no call to RNGkind or getRNG 
###% were done thereafter.
###% 
###% @return an object of class NativeSymbolInfo or NULL if no hook were found
###% 
RNGlib <- function(PACKAGE='', full=FALSE, hook="user_unif_rand", ...){
	
	if( !missing(PACKAGE) )
		full = TRUE
	if( !missing(hook) )
		hook <- match.arg(hook, c('user_unif_rand', 'user_unif_init', 'user_unif_nseed', 'user_unif_seedloc'))
	
	# lookup for the hook "user_unif_rand" in all the loaded libraries
	symb.unif_rand <- try( getNativeSymbolInfo(hook, PACKAGE=PACKAGE, ...), silent=TRUE)
	if( is(symb.unif_rand, 'try-error') ){
		
		if( !full ) '' else NULL
		
	}else if( PACKAGE=='' && is.null(symb.unif_rand$package) ){ 
		#special case for MS Windows when PACKAGE is not specified: if two 
		# RNGlibs are loaded, the first one is seen, not the last one as on Unix
		libs <- RNGlibs(full=TRUE, unlist=FALSE, hook=hook)
		w <- which(sapply(libs, function(l) identical(l$address, symb.unif_rand$address)))
		
		# returns full info or just the name
		if( full ) libs[[w]]
		else names(libs)[w]
		
	}else if( full ) symb.unif_rand
	else symb.unif_rand$package[['name']]
}

###% Returns the package that provides the current RNG managed by rstream
###% 
###% It returns the name of the package to which are currently passed the RNG 
###% calls (runif, set.seed).
###% This is either 'base' if core RNG is in use (e.g. Mersenne-Twister, Marsaglia-Multicarry, etc...) 
###% or the package that provides the actual RNG hooks called by the rstream 
###% wrapper hooks. This one was set either explicitly via RNGkind or implicitly 
###% when rstream was first loaded. In this latter case, the provider was identified 
###% at loading time as 'base' if core RNGs were in use or as the package that was 
###% providing the RNG hook 'user_unif_rand' if the RNG in used was "user-supplied".       
###%
RNGprovider <- function(kind=RNGkind(), user.supplied=FALSE){
	
	if( kind[1L] == 'user-supplied' || user.supplied ) RNGlib()		
	else 'base'
}

#' Directly Getting or Setting the RNG Seed
#' 
#' \code{RNGseed} directly gets/sets the current RNG seed \code{.Random.seed}.
#' It can typically be used to backup and restore the RNG state on exit of 
#' functions, enabling local RNG changes.
#' 
#' @param seed an RNG seed, i.e. an integer vector.
#' No validity check is performed, so it \strong{must} be a 
#' valid seed.
#' 
#' @return invisibly the current RNG seed when called with no arguments,
#' or the -- old -- value of the seed before changing it to 
#' \code{seed}. 
#' 
#' @export
#' @examples
#' 
#' # get current seed
#' RNGseed()
#' # directly set seed
#' old <- RNGseed(c(401L, 1L, 1L))
#' # show old/new seed description
#' showRNG(old)
#' showRNG()
#' 
#' # set bad seed
#' RNGseed(2:3)
#' try( showRNG() )
#' # recover from bad state
#' RNGrecovery()
#' showRNG()
#' 
#' # example of backup/restore of RNG in functions
#' f <- function(){
#' 	orng <- RNGseed()
#'  on.exit(RNGseed(orng))
#' 	RNGkind('Marsaglia')
#' 	runif(10)
#' }
#' 
#' sample(NA)
#' s <- .Random.seed
#' f()
#' identical(s, .Random.seed)
#' \dontshow{ stopifnot(identical(s, .Random.seed)) }
#' 
RNGseed <- function(seed){
	
	res <- if( missing(seed) ){
		if( exists('.Random.seed', where = .GlobalEnv) )
			get('.Random.seed', envir = .GlobalEnv)
	}else if( is.null(seed) ){
		if( exists('.Random.seed', where = .GlobalEnv) )
			rm('.Random.seed', envir = .GlobalEnv)
	}else{
		old <- RNGseed()
		assign('.Random.seed', seed, envir = .GlobalEnv)
		old
	}
	invisible(res)
}

#' \code{RNGrecovery} recovers from a broken state of \code{.Random.seed}, 
#' and reset the RNG settings to defaults.
#' 
#' @export
#' @rdname RNGseed
RNGrecovery <- function(){
	s <- as.integer(c(401,0,0))
	assign(".Random.seed", s, envir=.GlobalEnv)
	RNGkind("default", "default")
}

.getRNGattribute <- function(object){
	if( .hasSlot(object, 'rng') ) slot(object, 'rng')
	else if( .hasSlot(object, 'rng.seed') ) slot(object, 'rng.seed') # for back compatibility
	else if( !is.null(attr(object, 'rng')) ) attr(object, 'rng')
	else if( is.list(object) ){ # compatibility with package setRNG
		if( !is.null(object[['rng']]) ) object[['rng']]
		else if( is.list(object[['noise']]) && !is.null(object[['noise']][['rng']]) ) 
			object[['noise']][['rng']]
	}else NULL
}

#' Getting/Setting RNGs
#' 
#' \code{getRNG} returns the Random Number Generator (RNG) settings used for 
#' computing an object, using a suitable \code{.getRNG} S4 method to extract 
#' these settings.
#' For example, in the case of objects that result from multiple model fits, 
#' it would return the RNG settings used to compute the best fit.
#' 
#' This function handles single number RNG specifications in the following way:
#' \describe{
#' \item{integers}{Return them unchanged, considering them as encoded RNG kind 
#' specification (see \code{\link{RNG}}). No validity check is performed.}
#' \item{real numbers}{If \code{num.ok=TRUE} return them unchanged.
#' Otherwise, consider them as (pre-)seeds and pass them to \code{\link{set.seed}} 
#' to get a proper RNG seed.
#' Hence calling \code{getRNG(1234)} is equivalent to \code{set.seed(1234); getRNG()} 
#' (See examples).
#' }
#' }
#' 
#' @param object an R object from which RNG settings can be extracted, e.g. an 
#' integer vector containing a suitable value for \code{.Random.seed} or embedded 
#' RNG data, e.g., in S3/S4 slot \code{rng} or \code{rng$noise}.
#' @param ... extra arguments to allow extension and passed to a suitable S4 method 
#' \code{.getRNG} or \code{.setRNG}.
#' @param num.ok logical that indicates if single numeric (not integer) RNG data should be 
#' considered as a valid RNG seed (\code{TRUE}) or passed to \code{\link{set.seed}} 
#' into a proper RNG seed (\code{FALSE}) (See details and examples).
#' @param extract logical that indicates if embedded RNG data should be looked for and
#' extracted (\code{TRUE}) or if the object itself should be considered as an 
#' RNG specification.
#' @param recursive logical that indicates if embedded RNG data should be extracted 
#' recursively (\code{TRUE}) or only once (\code{FASE}). 
#' 
#' @return \code{getRNG}, \code{getRNG1}, \code{nextRNG} and \code{setRNG} 
#' usually return an integer vector of length > 2L, like \code{\link{.Random.seed}}.
#' 
#' \code{getRNG} and \code{getRNG1} return \code{NULL} if no RNG data was found.
#' 
#' @rdname rng
#' @seealso \code{\link{.Random.seed}}, \code{\link{showRNG}}
#' @export
#' 
#' @examples
#' # get current RNG settings
#' s <- getRNG()
#' head(s)
#' showRNG(s)
#' 
#' # get RNG from a given single numeric seed
#' s1234 <- getRNG(1234)
#' head(s1234)
#' showRNG(s1234)
#' # this is identical to the RNG seed as after set.seed()
#' set.seed(1234)
#' identical(s1234, .Random.seed)
#' # but if num.ok=TRUE the object is returned unchanged
#' getRNG(1234, num.ok=TRUE)
#' 
#' # single integer RNG data = encoded kind 
#' head(getRNG(1L))
#' 
#' # embedded RNG data
#' s <- getRNG(list(1L, rng=1234))
#' identical(s, s1234)
#'  
getRNG <- function(object, ..., num.ok=FALSE, extract=TRUE, recursive=TRUE){
	
	if( missing(object) || is.null(object) ) return( .getRNG() )
	
	# use RNG data from object if available
	if( extract && !is.null(rng <- .getRNGattribute(object)) ){
		if( recursive && hasRNG(rng) ) getRNG(rng, ..., num.ok=num.ok)
		else rng
	}else if( isNumber(object) ){
		if( num.ok && isReal(object) ) object
		else if( isInteger(object) ) object
		else nextRNG(object, ...) # return RNG as if after setting seed
	}else .getRNG(object, ...) # call S4 method on object
	
}

#' \code{hasRNG} tells if an object has embedded RNG data.
#' @export
#' @rdname rng
#' 
#' @examples
#' # test for embedded RNG data
#' hasRNG(1)
#' hasRNG( structure(1, rng=1:3) )
#' hasRNG( list(1, 2, 3) )
#' hasRNG( list(1, 2, 3, rng=1:3) )
#' hasRNG( list(1, 2, 3, noise=list(1:3, rng=1)) )
#' 
hasRNG <- function(object){
	!is.null(.getRNGattribute(object))
} 

#' \code{.getRNG} is an S4 generic that extract RNG settings from a variety of 
#' object types.
#' Its methods define the workhorse functions that are called by \code{getRNG}.
#' 
#' @rdname rng
#' @inline
#' @export
setGeneric('.getRNG', function(object, ...) standardGeneric('.getRNG') )
#' Default method that tries to extract RNG information from \code{object}, by 
#' looking sequentially to a slot named \code{'rng'}, a slot named \code{'rng.seed'}
#' or an attribute names \code{'rng'}.
#' 
#' It returns \code{NULL} if no RNG data was found.
setMethod('.getRNG', 'ANY',
	function(object, ...){
		.getRNGattribute(object)
	}
)
#' Returns the current RNG settings.
setMethod('.getRNG', 'missing',
	function(object){
		
		# return current value of .Random.seed
		# ensuring it exists first 
		if( !exists('.Random.seed', envir = .GlobalEnv) ) 
			sample(NA)
		
		return( get('.Random.seed', envir = .GlobalEnv) )
		
	}
)

#' Method for S3 objects, that aims at reproducing the behaviour of the function 
#' \code{getRNG} of the package \code{getRNG}. 
#' 
#' It sequentially looks for RNG data in elements \code{'rng'}, \code{noise$rng} 
#' if element \code{'noise'} exists and is a \code{list}, or in attribute \code{'rng'}.
#'  
setMethod('.getRNG', 'list',
	function(object){
		# lookup for some specific elements
		if( !is.null(object$rng) ) object$rng  
		else if( is.list(object$noise) ) object$noise$rng
		else attr(object, 'rng')
	}
)
#setMethod('.getRNG', 'rstream',
#		function(object){
#			object	
#		}
#)
#' Method for numeric vectors, which returns the object itself, coerced into an integer 
#' vector if necessary, as it is assumed to already represent a value for 
#' \code{\link{.Random.seed}}.
#' 
setMethod('.getRNG', 'numeric',
	function(object, ...){
		as.integer(object)
	}
)

#' \code{getRNG1} is an S4 generic that returns the \strong{initial} RNG settings 
#' used for computing an object.
#' For example, in the case of results from multiple model fies, it would 
#' return the RNG settings used to compute the \emph{first} fit.
#' 
#' \code{getRNG1} is defined to provide separate access to the RNG settings as 
#' they were at the very beginning of a whole computation, which might differ 
#' from the RNG settings returned by \code{getRNG}, that allows to reproduce the  
#' result only.
#' 
#' Think of a sequence of separate computations, from which only one result is 
#' used for the result (e.g. the one that maximises a likelihood): 
#' \code{getRNG1} would return the RNG settings to reproduce the complete sequence
#' of computations, while \code{getRNG} would return the RNG settings necessary to 
#' reproduce only the computation whose result has maximum likelihood.  
#' 
#' @rdname rng
#' @export
#' @inline
#' 
setGeneric('getRNG1', function(object, ...) standardGeneric('getRNG1') )
#' Default method that is identical to \code{getRNG(object, ...)}.
setMethod('getRNG1', 'ANY',
	function(object, ...){
		getRNG(object, ...)
	}
)

#' \code{nextRNG} returns the RNG settings as they would be after seeding with 
#' \code{object}.
#' 
#' @param ndraw number of draws to perform before returning the RNG seed.
#' 
#' @rdname rng
#' @export
#' @examples 
#' head(nextRNG())
#' head(nextRNG(1234))
#' head(nextRNG(1234, ndraw=10))
#' 
nextRNG <- function(object, ..., ndraw=0L){

	# get/restore .Random.seed on.exit
	orseed <- RNGseed()
	on.exit(RNGseed(orseed))
	
	# return next state of current RNG if object is missing
	if( missing(object) ){
		runif(1)
		return( getRNG() )
	}
	
	# extract RNG from object
	rng <- .getRNGattribute(object)
	if( !is.null(rng) ){
		on.exit()
		return( nextRNG(rng, ...) )
	}
	
	# only work for numeric seeds
	if( !is.numeric(object) )
		stop("Invalid seed: expecting a numeric seed.")
	
	# set RNG 
	.setRNG(object, ...)
	
	# perform some draws
	if( ndraw > 0 ){
		if( !isNumber(ndraw) )
			stop("Invalid value in argument `ndraw`: single numeric value expected.")
		runif(ndraw)
	}
	# return new RNG settings
	res <- RNGseed()
	res
}

.collapse <- function(x, sep=', ', n=7L){
	
	res <- paste(head(x, n), collapse=', ')
	if( length(x) > n )
		res <- paste(res, '...', sep=', ')
	res
}

#' \code{setRNG} set the current RNG with a seed, 
#' using a suitable \code{.setRNG} method to set these settings.
#'
#' @param check logical that indicates if only valid RNG kinds should be
#' accepted, or if invalid values should just throw a warning.
#' Note that this argument is used only on R >= 3.0.2.
#' 
#' @return \code{setRNG} invisibly returns the old RNG settings as 
#' they were before changing them.
#' 
#' @export 
#' @rdname rng
#' @examples 
#' 
#' obj <- list(x=1000, rng=123)
#' setRNG(obj)
#' rng <- getRNG()
#' runif(10)
#' set.seed(123)
#' rng.equal(rng)
#' 
setRNG <- function(object, ..., verbose=FALSE, check = TRUE){
	
	# do nothing if null
	if( is.null(object) ) return()
	
	# use RNG data from object if available
	rng <- getRNG(object, ...)
	if( !is.null(rng) && !identical(rng, object) ) return( setRNG(rng, ...) )
	
	# get/restore .Random.seed on.exit in case of errors
	orseed <- getRNG()
	on.exit({
		message("Restoring RNG settings probably due to an error in setRNG")
		RNGseed(orseed) 
	})

	# call S4 method on object
    # check validity of the seed
	tryCatch(.setRNG(object, ...)
            , warning = function(err){
                if( check && testRversion('> 3.0.1') 
                        && grepl("\\.Random\\.seed.* is not a valid", err$message) ){
                    stop("setRNG - Invalid RNG kind [", str_out(object), "]: "
							, err$message, '.'
							, call.=FALSE)
                }else{
                    warning(err)
                }
            } 
    )
	
	# cancel RNG restoration
	on.exit()
	if( verbose ) showRNG()			
	
	invisible(orseed)
}

#' \code{.setRNG} is an S4 generic that sets the current RNG settings, from a 
#' variety of specifications.
#' Its methods define the workhorse functions that are called by \code{setRNG}.
#' 
#' @inline
#' @rdname rng
#' @export 
setGeneric('.setRNG', function(object, ...) standardGeneric('.setRNG') )
#' Sets the RNG to kind \code{object}, assuming is a valid RNG kind:
#' it is equivalent to \code{RNGkind(object, ...}.
#' All arguments in \code{...} are passed to \code{\link{RNGkind}}.
#' 
#' @param verbose a logical that indicates if the new RNG settings should
#' be displayed.
#' 
#' @examples
#' # set RNG kind
#' old <- setRNG('Marsaglia')
#' # restore
#' setRNG(old)
setMethod('.setRNG', 'character',
	function(object, ...){
		if( length(object) == 1L )
			RNGkind(kind=object, ...)
		else
			RNGkind(kind=object[1L], normal.kind=object[2L])
	}
)

#' Sets the RNG settings using \code{object} directly the new value for 
#' \code{.Random.seed} or to initialise it with \code{\link{set.seed}}.
#' 
#' @examples 
#' 
#' # directly set .Random.seed
#' rng <- getRNG()
#' r <- runif(10)
#' setRNG(rng)
#' rng.equal(rng)
#' 
#' # initialise from a single number (<=> set.seed)
#' setRNG(123)
#' rng <- getRNG()
#' runif(10)
#' set.seed(123)
#' rng.equal(rng)
#' 
setMethod('.setRNG', 'numeric',
	function(object, ...){
		
		if( length(object) == 1L ){
			if( is.integer(object) ){ # set kind and draw once to fix seed
				RNGseed(object)
				# check validity of the seed
				tryCatch(runif(1)
					, error = function(err){					
						stop("setRNG - Invalid RNG kind [", object, "]: "
								, err$message, '.'
								, call.=FALSE)
				    }
                )
                RNGseed()
			}else{ # pass to set.seed
				set.seed(object, ...)
			}
		}else{
			seed <- as.integer(object)
			RNGseed(seed)
			# check validity of the seed
			tryCatch(runif(1)
			, error=function(err){					
				stop("setRNG - Invalid numeric seed ["
					, .collapse(seed, n=5), "]: ", err$message, '.'
					, call.=FALSE)
			    }
            )
			RNGseed(seed)			
		}
	}
)

#' \code{RNGdigest} computes a hash from the RNG settings associated with an 
#' object.
#' 
#' @import digest
#' @rdname RNGstr
#' @export
#' 
#' @examples 
#' # compute digest hash from RNG settings
#' RNGdigest()
#' RNGdigest(1234)
#' # no validity check is performed
#' RNGdigest(2:3)
#' 
RNGdigest <- function(object=getRNG()){
	
	x <- object
	object <- getRNG(x)
	
	# exit if no RNG was extracted
	if( is.null(object) ){
		warning("Found no embedded RNG data in object [", class(x),"]: returned NULL digest [", digest(NULL), '].')
		return(digest(NULL)) # TODO: return NULL
	}
		
	digest(object)
	
}

#' Comparing RNG Settings
#' 
#' \code{rng.equal} compares the RNG settings associated with two objects.
#' 
#' These functions return \code{TRUE} if the RNG settings are identical, 
#' and \code{FALSE} otherwise. 
#' The comparison is made between the hashes returned by \code{RNGdigest}.
#' 
#' @param x objects from which RNG settings are extracted
#' @param y object from which RNG settings are extracted
#' 
#' @return \code{rng.equal} and \code{rng.equal1} return a \code{TRUE} or 
#' \code{FALSE}.
#' 
#' @rdname rngcmp
#' @export
rng.equal <- function(x, y){
	if( missing(y) )
		y <- getRNG()
	identical(RNGdigest(x), RNGdigest(y))
}

#' \code{rng1.equal} tests whether two objects have identical 
#' \strong{initial} RNG settings.
#' 
#' @rdname rngcmp
#' @export
rng1.equal <- function(x, y){
	if( missing(y) )
		y <- getRNG()
	rng.equal(getRNG1(x), getRNG1(y))
}
