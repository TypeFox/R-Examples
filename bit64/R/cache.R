# /*
# R-Code for caching
# S3 atomic 64bit integers for R
# (c) 2011 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2011-12-11
# Last changed:  2011-12-11
# */

#! \name{cache}
#! \alias{cache}
#! \alias{newcache}
#! \alias{jamcache}
#! \alias{setcache}
#! \alias{getcache}
#! \alias{remcache}
#! \alias{print.cache}
#! \alias{still.identical}
#! \title{
#! 	Atomic Caching
#! }
#! \description{
#! 	Functions for caching results attached to atomic objects
#! }
#! \usage{
#! newcache(x)
#! jamcache(x)
#! cache(x)
#! setcache(x, which, value)
#! getcache(x, which)
#! remcache(x)
#! \method{print}{cache}(x, all.names = FALSE, pattern, \dots)
#! still.identical(x, y)
#! }
#! \arguments{
#!   \item{x}{
#!   an integer64 vector (or a cache object in case of \code{print.cache})
#! }
#!   \item{y}{
#!   an integer64 vector
#! }
#!   \item{which}{
#!   A character naming the object to be retrieved from the cache or to be stored in the cache
#! }
#!   \item{value}{
#!   An object to be stored in the cache 
#! }
#!   \item{all.names}{
#!   passed to \code{\link{ls}} when listing the cache content
#! }
#!   \item{pattern}{
#!   passed to \code{\link{ls}} when listing the cache content
#! }
#!   \item{\dots}{
#! 	ignored
#! }
#! }
#! \details{
#! 	A \code{cache} is an \code{link{environment}} attached to an atomic object with the \code{link{attrib}} name 'cache'. 
#! 	It contains at least a reference to the atomic object that carries the cache. 
#! 	This is used when accessing the cache to detect whether the object carrying the cache has been modified meanwhile.
#! 	Function \code{still.identical(x,y)} checks whether the objects \code{x} and \code{y} \cr
#! 	Function \code{newcache(x)} creates a new cache referencing  \code{x} \cr
#! 	Function \code{jamcache(x)} forces \code{x} to have a cache \cr
#! 	Function \code{cache(x)} returns the cache attached to \code{x} if it is not found to be outdated \cr
#! 	Function \code{setcache(x, which, value)} assigns a value into the cache of \code{x} \cr
#! 	Function \code{getcache(x, which)} gets cache value 'which' from \code{x} \cr
#! 	Function \code{remcache} removes the cache from \code{x} \cr
#! }
#! \value{
#! 	see details
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#! 	Functions that get and set small cache-content automatically when a cache is present: \code{\link{na.count}}, \code{\link{nvalid}}, \code{\link{is.sorted}}, \code{\link{nunique}} and \code{\link{nties}} \cr
#! 	Setting big caches with a relevant memory footprint requires a conscious decision of the user: \code{\link{hashcache}}, \code{\link{sortcache}}, \code{\link{ordercache}} and \code{\link{sortordercache}} \cr
#! 	Functions that use big caches: \code{\link{match.integer64}}, \code{\link{\%in\%.integer64}}, \code{\link{duplicated.integer64}}, \code{\link{unique.integer64}}, \code{\link{unipos}}, \code{\link{table.integer64}}, \code{\link{as.factor.integer64}}, \code{\link{as.ordered.integer64}}, \code{\link{keypos}}, \code{\link{tiepos}}, \code{\link{rank.integer64}}, \code{\link{prank}}, \code{\link{qtile}}, \code{\link{quantile.integer64}}, \code{\link{median.integer64}} and \code{\link{summary.integer64}} \cr
#! }
#! \examples{
#! 	x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! 	y <- x
#! 	still.identical(x,y)
#! 	y[1] <- NA
#! 	still.identical(x,y)
#! 	mycache <- newcache(x)
#! 	ls(mycache)
#! 	mycache
#! 	rm(mycache)
#! 	jamcache(x)
#! 	cache(x)
#! 	x[1] <- NA
#! 	cache(x)
#! 	getcache(x, "abc")
#! 	setcache(x, "abc", 1)
#! 	getcache(x, "abc")
#! 	remcache(x)
#! 	cache(x)
#! }
#! \keyword{ environment }

still.identical <- function(x, y){
  .Call("r_ram_truly_identical", x = x, y = y, PACKAGE = "bit64")
}

newcache <- function(x){
	env <- new.env()
	vmode <- typeof(x)
	if (vmode=="double" && is.integer64(x))
	  vmode <- "integer64"
	setattr(env, "class", c(paste("cache", vmode, sep="_"),"cache","environment"))
	assign("x", x, envir=env)
	env
}

jamcache <- function(x){
	cache <- attr(x, "cache")
	if (is.null(cache)){
		cache <- newcache(x)
		setattr(x, "cache", cache)
	}else
		if (!still.identical(x, get("x", envir=cache, inherits=FALSE))){
			cache <- newcache(x)
			setattr(x, "cache", cache)
			warning("replaced outdated cache with empty cache")
		}
	cache
}

cache <- function(x){
	cache <- attr(x, "cache")
	if (is.null(cache) || still.identical(x, get("x", envir=cache, inherits=FALSE)))
		cache
	else{ 
		remcache(x)
		warning("removed outdated cache")
		NULL
	}
}

setcache <- function(x, which, value){
	  env <- jamcache(x)
	  assign(which, value, envir=env)
	  env
}

getcache <- function(x, which){
	cache <- attr(x, "cache")
	if (is.null(cache))
	  return(NULL)
	if (still.identical(x, get("x", envir=cache, inherits=FALSE))){
		if (exists(which, envir=cache, inherits=FALSE))
			get(which, envir=cache, inherits=FALSE)
		else
			NULL
	}else{ 
		remcache(x)
		warning("removed outdated cache")
		NULL
	}
}

remcache <- function(x){
		setattr(x, "cache", NULL)
	invisible()
}

print.cache<- function(x, all.names=FALSE, pattern, ...){
  l <- ls(x, all.names, pattern=pattern)
  cat(class(x)[1], ": ", paste(l, collapse=" - "), "\n", sep="")
  invisible(l)
}


#! \name{hashcache}
#! \alias{hashcache}
#! \alias{sortcache}
#! \alias{sortordercache}
#! \alias{ordercache}
#! \title{
#! 		Big caching of hashing, sorting, ordering
#! }
#! \description{
#! 	Functions to create cache that accelerates many operations
#! }
#! \usage{
#! hashcache(x, nunique=NULL, \dots)
#! sortcache(x, has.na = NULL)
#! sortordercache(x, has.na = NULL, stable = NULL)
#! ordercache(x, has.na = NULL, stable = NULL, optimize = "time")
#! }
#! \arguments{
#!   \item{x}{
#! 		an atomic vector (note that currently only integer64 is supported)
#! }
#!   \item{nunique}{ giving \emph{correct} number of unique elements can help reducing the size of the hashmap }
#!   \item{has.na}{
#! boolean scalar defining whether the input vector might contain \code{NA}s. If we know we don't have NAs, this may speed-up.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{stable}{
#! boolean scalar defining whether stable sorting is needed. Allowing non-stable may speed-up.
#! }
#!   \item{optimize}{
#! by default ramsort optimizes for 'time' which requires more RAM,
#! set to 'memory' to minimize RAM requirements and sacrifice speed
#! }
#!   \item{\dots}{
#! 		passed to \code{\link{hashmap}}
#! }
#! }
#! \details{
#! 	The result of relative expensive operations \code{\link{hashmap}}, \code{\link{ramsort}}, \code{\link{ramsortorder}} and \code{\link{ramorder}} can be stored in a cache in order to avoid multiple excutions. Unless in very specific situations, the recommended method is \code{hashsortorder} only.
#! }
#! \note{
#!   Note that we consider storing the big results from sorting and/or ordering as a relevant side-effect, 
#! and therefore storing them in the cache should require a conscious decision of the user.
#! }
#! \value{
#! 	\code{x} with a \code{\link{cache}} that contains the result of the expensive operations, possible together with small derived information (such as \code{\link{nunique.integer64}}) and previously cached results.
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#! 	\code{\link{cache}} for caching functions and \code{\link{nunique}} for methods bennefitting from small caches
#! }
#! \examples{
#! 	x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#!  sortordercache(x)
#! }
#! \keyword{ environment }

hashcache <-function(x, nunique=NULL, ...){
	env <- jamcache(x)
	if (is.null(nunique))
		nunique <- env$nunique
	env <- hashmap(x, nunique=nunique, cache=env, ...)
	if (is.null(nunique) && env$nunique<sqrt(length(x)))
		env <- hashmap(x, nunique=env$nunique, cache=env, ...)
	na.count(x) # since x has cache, na.count() will update the cache, unless its already there
	# different from sortcache, ordercache and sortordercache we do not set nties: hastab is too expensive
	invisible(env)
}

sortcache <- function(x, has.na = NULL){
	if (is.null(has.na)){
		na.count <- getcache(x, "na.count")
		if (is.null(na.count))
			has.na <- TRUE
		else
			has.na <- na.count > 0
	}
	s <- clone(x)
    na.count <- ramsort(s, has.na = has.na, na.last = FALSE, decreasing = FALSE, stable = FALSE, optimize = "time")
	nut <- .Call("r_ram_integer64_sortnut", x = s, PACKAGE = "bit64")
    setcache(x, "sort", s)
    setcache(x, "na.count", na.count)
    setcache(x, "nunique", nut[[1]])
    setcache(x, "nties", nut[[2]])
	invisible(x)
}


sortordercache <- function(x, has.na = NULL, stable = NULL){
	if (is.null(has.na)){
		na.count <- getcache(x, "na.count")
		if (is.null(na.count))
			has.na <- TRUE
		else
			has.na <- na.count > 0
	}
	if (is.null(stable)){
		nunique <- getcache(x, "nunique")
		if (is.null(nunique))
		  stable <- TRUE
		else
		  stable <- nunique < length(x)
	}
	s <- clone(x)
	o <- seq_along(x)
    na.count <- ramsortorder(s, o, has.na = has.na, na.last = FALSE, decreasing = FALSE, stable = stable, optimize = "time")
	nut <- .Call("r_ram_integer64_sortnut", x = s, PACKAGE = "bit64")
    setcache(x, "sort", s)
    setcache(x, "order", o)
    setcache(x, "na.count", na.count)
    setcache(x, "nunique", nut[[1]])
    setcache(x, "nties", nut[[2]])
	invisible(x)
}


ordercache <- function(x, has.na = NULL, stable = NULL, optimize = "time"){
	if (is.null(has.na)){
		na.count <- getcache(x, "na.count")
		if (is.null(na.count))
			has.na <- TRUE
		else
			has.na <- na.count > 0
	}
	if (is.null(stable)){
		nunique <- getcache(x, "nunique")
		if (is.null(nunique))
		  stable <- TRUE
		else
		  stable <- nunique < length(x)
	}
	o <- seq_along(x)
    na.count <- ramorder(x, o, has.na = has.na, na.last = FALSE, decreasing = FALSE, stable = stable, optimize = optimize)
	nut <- .Call("r_ram_integer64_ordernut", table = x, order = o, PACKAGE = "bit64")
    setcache(x, "order", o)
    setcache(x, "na.count", na.count)
    setcache(x, "nunique", nut[[1]])
    setcache(x, "nties", nut[[2]])
	invisible(x)
}



#! \name{is.sorted.integer64}
#! \alias{is.sorted.integer64}
#! \alias{na.count.integer64}
#! \alias{nvalid.integer64}
#! \alias{nunique.integer64}
#! \alias{nties.integer64}
#! \title{
#! 	Small cache access methods
#! }
#! \description{
#! 	These methods are packaged here for methods in packages \code{bit64} and \code{ff}.
#! }
#! \usage{
#! 	\method{is.sorted}{integer64}(x, \dots)
#! 	\method{na.count}{integer64}(x, \dots)
#! 	\method{nvalid}{integer64}(x, \dots)
#! 	\method{nunique}{integer64}(x, \dots)
#! 	\method{nties}{integer64}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{
#! 	some object
#! 	}
#!   \item{\dots}{
#! 	ignored
#! 	}
#! }
#! \details{
#!   All these functions benefit from a \code{\link{sortcache}}, \code{\link{ordercache}} or \code{\link{sortordercache}}.  
#!   \code{na.count}, \code{nvalid} and \code{nunique} also benefit from a \code{\link{hashcache}}.
#!	\cr
#! 	\code{is.sorted} checks for sortedness of \code{x} (NAs sorted first) \cr
#!  \code{na.count} returns the number of \code{NA}s \cr 
#!  \code{nvalid} returns the number of valid data points, usually \code{\link{length}} minus \code{na.count}. \cr
#!  \code{nunique} returns the number of unique values \cr
#!  \code{nties} returns the number of tied values. 
#! }
#! \note{
#! 	If a \code{\link{cache}} exists but the desired value is not cached, 
#!  then these functions will store their result in the cache. 
#!  We do not consider this a relevant side-effect, 
#!  since these small cache results do not have a relevant memory footprint.
#! }
#! \value{
#! 	\code{is.sorted} returns a logical scalar, the other methods return an integer scalar.
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#! 	\code{\link{cache}} for caching functions and \code{\link{sortordercache}} for functions creating big caches
#! }
#! \examples{
#! 	x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#!  length(x)
#!  na.count(x)
#!  nvalid(x)
#!  nunique(x)
#!  nties(x)
#!  table.integer64(x)
#!  x
#! }
#! \keyword{ environment }
#! \keyword{ methods }


na.count.integer64 <- function(x, ...){
  env <- cache(x)
  if (is.null(env)){
	.Call("r_ram_integer64_nacount", x = x, PACKAGE = "bit64")
  }else{
    if (exists("na.count", envir=env, inherits=FALSE))
		get("na.count", envir=env, inherits=FALSE)
	else{
		ret <- .Call("r_ram_integer64_nacount", x = x, PACKAGE = "bit64")
		assign("na.count", ret, envir=env)
		ret
	}
  }
}

nvalid.integer64 <- function(x, ...){
	length(x) - na.count(x)
}

is.sorted.integer64 <- function(x, ...){
  env <- cache(x)
  if (is.null(env)){
	.Call("r_ram_integer64_issorted_asc", x = x, PACKAGE = "bit64")
  }else{
    if (exists("is.sorted", envir=env, inherits=FALSE))
		get("is.sorted", envir=env, inherits=FALSE)
	else{
		ret <- .Call("r_ram_integer64_issorted_asc", x = x, PACKAGE = "bit64")
		assign("is.sorted", ret, envir=env)
		ret
	}
  }
}


nunique.integer64 <- function(x, ...){
	env <- cache(x)
	if(is.null(env))
		has.cache <- FALSE
	else{
		if (exists("nunique", envir=env, inherits=FALSE))
			return(get("nunique", envir=env, inherits=FALSE))
		else
			has.cache <- TRUE
	}
	if (is.sorted(x)){
		ret <- .Call("r_ram_integer64_sortnut"
		, x = x
		, PACKAGE = "bit64"
		)
		if (has.cache){
			assign("nunique", ret[1], envir=env)
			assign("nties", ret[2], envir=env)
		}
		ret[1]
	}else{
		h <- hashmap(x)
		if (has.cache)
		  assign("nunique", h$nunique, envir=env)
		h$nunique
	}
}

nties.integer64 <- function(x, ...){
	cv <- getcache(x, "nties")
	if (is.null(cv)){
		if (is.sorted(x)){
			cv <- .Call("r_ram_integer64_sortnut"
			, x = x
			, PACKAGE = "bit64"
			)[2]
		}else{
		    s <- clone(x)
			na.count <- ramsort(s, has.na = TRUE, na.last = FALSE, decreasing = FALSE, stable = FALSE, optimize = "time")
			cv <- .Call("r_ram_integer64_sortnut", x = s, PACKAGE = "bit64")[[2]]
		}
	}
	cv
}

