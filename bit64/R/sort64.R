# /*
# R-Code for sorting and ordering
# S3 atomic 64bit integers for R
# (c) 2011 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2011-12-11
# Last changed:  2011-12-11
# */

#! \name{ramsort.integer64}
#! \alias{ramsort.integer64}
#! \alias{shellsort.integer64}
#! \alias{quicksort.integer64}
#! \alias{mergesort.integer64}
#! \alias{radixsort.integer64}
#! \alias{ramorder.integer64}
#! \alias{shellorder.integer64}
#! \alias{quickorder.integer64}
#! \alias{mergeorder.integer64}
#! \alias{radixorder.integer64}
#! \alias{ramsortorder.integer64}
#! \alias{shellsortorder.integer64}
#! \alias{quicksortorder.integer64}
#! \alias{mergesortorder.integer64}
#! \alias{radixsortorder.integer64}
#! \title{
#!    Low-level intger64 methods for in-RAM sorting and ordering
#! }
#! \description{
#!   Fast low-level methods for sorting and ordering. 
#!   The \code{..sortorder} methods do sorting and ordering at once, which requires more RAM than ordering but is (almost) as fast as as sorting.
#! }
#! \note{
#!  Note that these methods purposely violate the functional programming paradigm: they are called for the side-effect of changing some of their arguments.
#!  The \code{sort}-methods change \code{x}, the \code{order}-methods change \code{i}, and the \code{sortoder}-methods change both \code{x} and \code{i}
#! }
#! \usage{
#! \method{shellsort}{integer64}(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE, \dots)
#! \method{shellsortorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, \dots)
#! \method{shellorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, \dots)
#! \method{mergesort}{integer64}(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE, \dots)
#! \method{mergeorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, \dots)
#! \method{mergesortorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, \dots)
#! \method{quicksort}{integer64}(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE
#! , restlevel=floor(1.5*log2(length(x))), \dots)
#! \method{quicksortorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE
#! , restlevel=floor(1.5*log2(length(x))), \dots)
#! \method{quickorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE
#! , restlevel=floor(1.5*log2(length(x))), \dots)
#! \method{radixsort}{integer64}(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE, radixbits=8L, \dots)
#! \method{radixsortorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, radixbits=8L, \dots)
#! \method{radixorder}{integer64}(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, radixbits=8L, \dots)
#! \method{ramsort}{integer64}(x, has.na = TRUE, na.last=FALSE, decreasing = FALSE, stable = TRUE
#! , optimize = c("time", "memory"), VERBOSE = FALSE, \dots)
#! \method{ramsortorder}{integer64}(x, i, has.na = TRUE, na.last=FALSE, decreasing = FALSE, stable = TRUE
#! , optimize = c("time", "memory"), VERBOSE = FALSE, \dots)
#! \method{ramorder}{integer64}(x, i, has.na = TRUE, na.last=FALSE, decreasing = FALSE, stable = TRUE
#! , optimize = c("time", "memory"), VERBOSE = FALSE, \dots)
#! }
#! \arguments{
#!   \item{x}{ a vector to be sorted by \code{\link{ramsort}} and \code{\link{ramsortorder}}, i.e. the output of  \code{\link{sort}} }
#!   \item{i}{ integer positions to be modified by \code{\link{ramorder}} and \code{\link{ramsortorder}}, default is 1:n, in this case the output is similar to \code{\link{order}} }
#!   \item{has.na}{
#! boolean scalar defining whether the input vector might contain \code{NA}s. If we know we don't have NAs, this may speed-up.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{na.last}{
#! boolean scalar telling ramsort whether to sort \code{NA}s last or first.
#! \emph{Note} that 'boolean' means that there is no third option \code{NA} as in \code{\link{sort}}
#! }
#!   \item{decreasing}{
#! boolean scalar telling ramsort whether to sort increasing or decreasing
#! }
#!   \item{stable}{
#! boolean scalar defining whether stable sorting is needed. Allowing non-stable may speed-up.
#! }
#!   \item{optimize}{
#! by default ramsort optimizes for 'time' which requires more RAM,
#! set to 'memory' to minimize RAM requirements and sacrifice speed
#! }
#!   \item{restlevel}{
#! number of remaining recursionlevels before \code{quicksort} switches from recursing to \code{shellsort}
#! }
#!   \item{radixbits}{
#! 	size of radix in bits
#! }
#!   \item{VERBOSE}{
#!   cat some info about chosen method
#! }
#!   \item{\dots}{ further arguments, passed from generics, ignored in methods }
#! }
#! \details{
#!  see \code{\link[bit]{ramsort}}
#! }
#! \value{
#!   These functions return the number of \code{NAs} found or assumed during sorting
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ programming }
#! \keyword{ manip }
#! \seealso{ \code{\link{ramsort}} for the generic, \code{\link[ff]{ramsort.default}} for the methods provided by package \code{\link[ff]{ff}}, \code{\link{sort.integer64}} for the sort interface and \code{\link{sortcache}} for caching the work of sorting}
#! \examples{
#!   x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#!   x
#!   message("ramsort example")
#!   s <- clone(x)
#!   ramsort(s)
#!   message("s has been changed in-place - whether or not ramsort uses an in-place algorithm")
#!   s
#!   message("ramorder example")
#!   s <- clone(x)
#!   o <- seq_along(s)
#!   ramorder(s, o)
#!   message("o has been changed in-place - s remains unchanged")
#!   s
#!   o
#!   s[o]
#!   message("ramsortorder example")
#!   o <- seq_along(s)
#!   ramsortorder(s, o)
#!   message("s and o have both been changed in-place - this is much faster")
#!   s
#!   o
#! }

shellsort.integer64 <- function(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE, ...)
{
  force(x)
  .Call("r_ram_integer64_shellsort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "bit64"
  )
}
shellsortorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, ...)
{
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  .Call("r_ram_integer64_shellsortorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "bit64"
  )
}
shellorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, ...)
{
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  .Call("r_ram_integer64_shellorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "bit64"
  )
}

mergesort.integer64 <- function(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE, ...){
  force(x)
  .Call("r_ram_integer64_mergesort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "bit64"
  )
}

mergeorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, ...){
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  .Call("r_ram_integer64_mergeorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "bit64"
  )
}

mergesortorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE, ...){
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  .Call("r_ram_integer64_mergesortorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "bit64"
  )
}


quicksort.integer64 <- function(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE
, restlevel=floor(1.5*log2(length(x)))
, ...){
  force(x)
  if (restlevel<0)
    restlevel = 0L
  .Call("r_ram_integer64_quicksort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , restlevel = as.integer(restlevel)
  , PACKAGE = "bit64"
  )
}

quicksortorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE
, restlevel=floor(1.5*log2(length(x)))
, ...){
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  if (restlevel<0)
    restlevel = 0L
  .Call("r_ram_integer64_quicksortorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , restlevel = as.integer(restlevel)
  , PACKAGE = "bit64"
  )
}

quickorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE
, restlevel=floor(1.5*log2(length(x)))
, ...){
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  if (restlevel<0)
    restlevel = 0L
  .Call("r_ram_integer64_quickorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , restlevel = as.integer(restlevel)
  , PACKAGE = "bit64"
  )
}

radixsort.integer64 <- function(x, has.na=TRUE, na.last=FALSE, decreasing=FALSE
, radixbits=8L
, ...)
{
  stopifnot(radixbits %in% c(1L, 2L, 4L, 8L, 16L))
  force(x)
  .Call("r_ram_integer64_radixsort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , radixbits = as.integer(radixbits)
  , PACKAGE = "bit64"
  )
}

radixsortorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE
, radixbits=8L
, ...)
{
  stopifnot(radixbits %in% c(1L, 2L, 4L, 8L, 16L))
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  .Call("r_ram_integer64_radixsortorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , radixbits = as.integer(radixbits)
  , PACKAGE = "bit64"
  )
}

radixorder.integer64 <- function(x, i, has.na=TRUE, na.last=FALSE, decreasing=FALSE
, radixbits=8L
, ...)
{
  stopifnot(radixbits %in% c(1L, 2L, 4L, 8L, 16L))
  force(x)
  force(i)
  if (!is.integer(i)) 
    stop("i must be integer")
  if (length(i) != length(x)) 
    stop("lengths of x and i don't match")  
  .Call("r_ram_integer64_radixorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , radixbits = as.integer(radixbits)
  , PACKAGE = "bit64"
  )
}

ramsort.integer64 <- function (x
, has.na = TRUE
, na.last=FALSE
, decreasing = FALSE
, stable = TRUE
, optimize = c("time", "memory")
, VERBOSE = FALSE
, ...
)
{
	optimize <- match.arg(optimize)
	if (is.null(names(x))){
		if (optimize == "time"){
			if (length(x)<2048){
				if (VERBOSE) 
					cat("ramsort selected mergesort\n")
				mergesort(x, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}else if (length(x)<16777216){
				if (VERBOSE) 
					cat("ramsort selected radix8sort\n")
				radixsort(x, radixbits=8L, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}else{
				if (VERBOSE) 
					cat("ramsort selected radix4sort\n")
				radixsort(x, radixbits=4L, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}
		}else{
			if (VERBOSE) 
				cat("ramsort selected quicksort\n")
			quicksort(x, has.na = has.na, na.last = na.last, decreasing = decreasing)
		}
	}else{
		if (stable || optimize == "time"){
			i <- seq_along(x)
		    if (length(x)<2048){
				if (VERBOSE) 
					cat("ramsortorder selected mergesortorder\n")
				ret <- mergesortorder(x, i, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}else if (length(x)<2097152){
				if (VERBOSE) 
					cat("ramsortorder selected radix8sortorder\n")
				ret <- radixsortorder(x, i, radixbits=8L, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}else{
				if (VERBOSE) 
					cat("ramsortorder selected radix4sortorder\n")
				ret <- radixsortorder(x, i, radixbits=4L, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}
		}else{
			if (VERBOSE) 
				cat("ramsort selected quicksortorder\n")
			i <- seq_along(x)
			ret <- quicksortorder(x, i, has.na = has.na, na.last = na.last, decreasing = decreasing)
		}
    setattr(x, "names", names(x)[i])
    ret
	}
}

ramsortorder.integer64 <- function (x
, i
, has.na = TRUE
, na.last=FALSE
, decreasing = FALSE
, stable = TRUE
, optimize = c("time", "memory")
, VERBOSE = FALSE
, ...
)
{
	optimize <- match.arg(optimize)
	if (is.null(names(x)) & is.null(names(i))){
		if (stable || optimize == "time") {
		    if (length(x)<2048){
				if (VERBOSE) 
					cat("ramsortorder selected mergesortorder\n")
				mergesortorder(x, i, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}else if (length(x)<16777216){
				if (VERBOSE) 
					cat("ramsortorder selected radix8sortorder\n")
				radixsortorder(x, i, radixbits=8L, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}else{
				if (VERBOSE) 
					cat("ramsortorder selected radix4sortorder\n")
				radixsortorder(x, i, radixbits=4L, has.na = has.na, na.last = na.last, decreasing = decreasing)
			}
		}else{
			if (VERBOSE) 
				cat("ramsortorder selected quicksortorder\n")
			quicksortorder(x, i, has.na = has.na, na.last = na.last, decreasing = decreasing)
		}
	}else
	  stop("names not supported")
}

ramorder.integer64 <- function (x
, i
, has.na = TRUE
, na.last=FALSE
, decreasing = FALSE
, stable = TRUE
, optimize = c("time", "memory")
, VERBOSE = FALSE
, ...
)
{
	optimize <- match.arg(optimize)
	if (is.null(names(x)) & is.null(names(i))){
		if (stable) {
			if (VERBOSE) 
				cat("ramorder selected mergeorder\n")
			mergeorder(x, i, has.na = has.na, na.last = na.last, decreasing = decreasing)
		}else{
			if (VERBOSE) 
				cat("ramorder selected quickorder\n")
			quickorder(x, i, has.na = has.na, na.last = na.last, decreasing = decreasing)
		}
	}else
	  stop("names not supported")
}


#! \name{sort.integer64}
#! \alias{sort.integer64}
#! \alias{order.integer64}
#! \title{
#!    High-level intger64 methods for sorting and ordering
#! }
#! \description{
#!   Fast high-level methods for sorting and ordering. 
#!   These are wrappers to \code{\link{ramsort}} and friends and do not modify their arguments.
#! }
#! \usage{
#! \method{sort}{integer64}(x, decreasing = FALSE, has.na = TRUE, na.last = TRUE, stable = TRUE
#! , optimize = c("time", "memory"), VERBOSE = FALSE, \dots)
#! \method{order}{integer64}(\dots, na.last = TRUE, decreasing = FALSE, has.na = TRUE, stable = TRUE
#! , optimize = c("time", "memory"), VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{x}{ a vector to be sorted by \code{\link{ramsort}} and \code{\link{ramsortorder}}, i.e. the output of  \code{\link{sort}} }
#!   \item{has.na}{
#! boolean scalar defining whether the input vector might contain \code{NA}s. If we know we don't have NAs, this may speed-up.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{na.last}{
#! boolean scalar telling ramsort whether to sort \code{NA}s last or first.
#! \emph{Note} that 'boolean' means that there is no third option \code{NA} as in \code{\link{sort}}
#! }
#!   \item{decreasing}{
#! boolean scalar telling ramsort whether to sort increasing or decreasing
#! }
#!   \item{stable}{
#! boolean scalar defining whether stable sorting is needed. Allowing non-stable may speed-up.
#! }
#!   \item{optimize}{
#! by default ramsort optimizes for 'time' which requires more RAM,
#! set to 'memory' to minimize RAM requirements and sacrifice speed
#! }
#!   \item{VERBOSE}{
#!   cat some info about chosen method
#! }
#!   \item{\dots}{ further arguments, passed from generics, ignored in methods }
#! }
#! \details{
#!  see \code{\link{sort}} and \code{\link{order}}
#! }
#! \value{
#!   \code{sort} returns the sorted vector and \code{vector} returns the order positions. 
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ programming }
#! \keyword{ manip }
#! \seealso{ \code{\link[=sort.integer64]{sort}}, \code{\link{sortcache}} }
#! \examples{
#!   x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#!   x
#!   sort(x)
#!   message("the following has default optimize='time' which is faster but requires more RAM
#! , this calls 'ramorder'")
#!   order.integer64(x)
#!   message("slower with less RAM, this calls 'ramsortorder'")
#!   order.integer64(x, optimize="memory")
#! }

if (FALSE){
	library(bit64)
	x <- as.integer64(c(sample(1e7),NA))
	#system.time(sortcache(x))[3]
	# system.time(ordercache(x))[3]
	system.time(sortordercache(x))[3]

	# system.time(s <- sort(x, na.last=FALSE, decreasing=FALSE))[3]
	# stopifnot(identical(s, {xs<-clone(x);ramsort(xs, na.last=FALSE, decreasing=FALSE);xs}))
	# system.time(s <- sort(x, na.last=TRUE, decreasing=FALSE))[3]
	# stopifnot(identical(s, {xs<-clone(x);ramsort(xs, na.last=TRUE, decreasing=FALSE);xs}))
	# system.time(s <- sort(x, na.last=FALSE, decreasing=TRUE))[3]
	# stopifnot(identical(s, {xs<-clone(x);ramsort(xs, na.last=FALSE, decreasing=TRUE);xs}))
	# system.time(s <- sort(x, na.last=TRUE, decreasing=TRUE))[3]
	# stopifnot(identical(s, {xs<-clone(x);ramsort(xs, na.last=TRUE, decreasing=TRUE);xs}))
	
	system.time(o <- order.integer64(x, na.last=FALSE, decreasing=FALSE))[3]
	stopifnot(identical(o, {xo<-seq_along(x);ramorder(x, xo, na.last=FALSE, decreasing=FALSE);xo}))
	system.time(o <- order.integer64(x, na.last=TRUE, decreasing=FALSE))[3]
	stopifnot(identical(o, {xo<-seq_along(x);ramorder(x, xo, na.last=TRUE, decreasing=FALSE);xo}))
	system.time(o <- order.integer64(x, na.last=FALSE, decreasing=TRUE))[3]
	stopifnot(identical(o, {xo<-seq_along(x);ramorder(x, xo, na.last=FALSE, decreasing=TRUE);xo}))
	system.time(o <- order.integer64(x, na.last=TRUE, decreasing=TRUE))[3]
	stopifnot(identical(o, {xo<-seq_along(x);ramorder(x, xo, na.last=TRUE, decreasing=TRUE);xo}))
	
}

sort.integer64 <- function(x
, decreasing = FALSE
, has.na = TRUE
, na.last = TRUE
, stable = TRUE
, optimize = c("time", "memory")
, VERBOSE = FALSE
, ...
){
  do.na.last <- is.na(na.last) || na.last
  c <- cache(x)
  if (!is.null(c$sort)){
		if (do.na.last || decreasing){
			s <- double(length(x))
			.Call("r_ram_integer64_sortsrt"
			, x = c$sort
			, na_count   = as.integer(na.count <- c$na.count)
			, na_last    = as.logical(do.na.last)
			, decreasing = as.logical(decreasing)
			, s		 	 = s
			, PACKAGE = "bit64"
			)
			setattr(s, "class", "integer64")
		}else 
			s <- c$sort  # here we save copying at all
  }else if (!is.null(c$order)){
		if (do.na.last || decreasing){
			s <- double(length(x))
			.Call("r_ram_integer64_sortsrt"
			, x = x[c$order]
			, na_count   = as.integer(na.count <- c$na.count)
			, na_last    = as.logical(do.na.last)
			, decreasing = as.logical(decreasing)
			, s		 	 = s
			, PACKAGE = "bit64"
			)
			setattr(s, "class", "integer64")
		}else 
			s <- x[c$order]
  }else{
    if (identical(c$na.count, 0L))
	  has.na <- FALSE
		s <- clone(x)
		na.count <- ramsort(
			s
		, has.na=has.na
		, na.last=do.na.last
		, decreasing=decreasing
		, stable=stable
		, optimize = optimize
		, VERBOSE = FALSE
		)
  }
  if (is.na(na.last) && na.count)
		length(s) <- length(s) - na.count
  s
}


order.integer64 <- function(
  ...
, na.last = TRUE
, decreasing = FALSE
, has.na = TRUE
, stable = TRUE
, optimize = c("time", "memory")
, VERBOSE = FALSE
){
  do.na.last <- is.na(na.last) || na.last
	# COPY ON MODIFY is broken for reading from list(...)
	# because list(...) creates a copy of all ... and this invalidates our caches
	# therefore we go this sick workaround
	argsymbols <- as.list(substitute(list(...)))[-1L]
	argframe <- parent.frame()
	A <- function(i)eval(argsymbols[[i]], argframe)
	N <- length(argsymbols)
  if (N!=1L)
	stop("can only order one vector at the moment")
  x <- A(1)
  c <- cache(x)
  if (!is.null(c$order)){
		if (do.na.last || decreasing){
			o <- integer(length(x))
			if (is.null(c$sort)){
				.Call("r_ram_integer64_orderord"
				, x = x
				, i = c$order
				, na_count   = as.integer(na.count <- c$na.count)
				, na_last    = as.logical(do.na.last)
				, decreasing = as.logical(decreasing)
				, o		 	 = o
				, PACKAGE = "bit64"
				)
			}else{
				.Call("r_ram_integer64_sortorderord"
				, x = c$sort
				, i = c$order
				, na_count   = as.integer(na.count <- c$na.count)
				, na_last    = as.logical(do.na.last)
				, decreasing = as.logical(decreasing)
				, o		 	 = o
				, PACKAGE = "bit64"
				)
			}
  		}else 
			o <- c$order  # here we save copying at all
  }else{
	  if (identical(c$na.count, 0L))
		has.na <- FALSE
	  optimize <- match.arg(optimize)
	  o <- seq_along(x)
	  if (optimize=="time"){
		  s <- clone(x)
		  na.count <- ramsortorder(s, o
		  , has.na=has.na
		  , na.last=do.na.last
		  , decreasing=decreasing
		  , stable=stable
		  , optimize = optimize
		  , VERBOSE = FALSE
		  )
	  }else{
		  na.count <- ramorder(x, o
		  , has.na=has.na
		  , na.last=do.na.last
		  , decreasing=decreasing
		  , stable=stable
		  , optimize = optimize
		  , VERBOSE = FALSE
		  )
	  }
	}
	if (is.na(na.last) && na.count)
	  length(o) <- length(o) - na.count
	o
}


