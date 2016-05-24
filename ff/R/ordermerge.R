# sorting module for ff
# (c) 2010 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2010-01-01
# Last changed: 2010-06-06

#source("c:/MWP/eAnalysis/ff/R/ordermerge.R")

#! \name{ramsort.default}
#! \alias{ramsort.default}
#! \alias{mergesort.default}
#! \alias{radixsort.default}
#! \alias{keysort.default}
#! \alias{shellsort.default}
#! \title{
#! Sorting: Sort R vector in-RAM and in-place
#! }
#! \description{
#! Function \code{ramsort} will sort the input vector in-place (without making a copy) and return the number of NAs found
#! }
#! \usage{
#! \method{ramsort}{default}(x, has.na = TRUE, na.last = TRUE, decreasing = FALSE
#! , optimize = c("time", "memory"), VERBOSE = FALSE, \dots)
#! \method{mergesort}{default}(x, has.na = TRUE, na.last = TRUE, decreasing = FALSE, \dots)
#! \method{radixsort}{default}(x, has.na = TRUE, na.last = TRUE, decreasing = FALSE, \dots)
#! \method{keysort}{default}(x, keyrange=range(x, na.rm=has.na), has.na = TRUE
#! , na.last = TRUE, decreasing = FALSE, \dots)
#! \method{shellsort}{default}(x, has.na = TRUE, na.last = TRUE, decreasing = FALSE, \dots)
#! }
#! \arguments{
#!   \item{x}{
#! an atomic R vector
#! }
#!   \item{keyrange}{
#! an integer vector with two values giving the smallest and largest possible value in x, note that you should give this explicitely for best performance, relying on the default needs one pass over the data to determine the range
#! }
#!   \item{has.na}{
#! boolean scalar telling ramsort whether the vector might contain \code{NA}s.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{na.last}{
#! boolean scalar telling ramsort whether to sort \code{NA}s last or first.
#! \emph{Note} that 'boolean' means that there is no third option \code{NA} as in \code{\link{sort}}
#! }
#!   \item{decreasing}{
#! boolean scalar telling ramsort whether to sort increasing or decreasing
#! }
#!   \item{optimize}{
#! by default ramsort optimizes for 'time' which requires more RAM,
#! set to 'memory' to minimize RAM requirements and sacrifice speed
#! }
#!   \item{VERBOSE}{
#!   cat some info about chosen method
#! }
#!   \item{\dots}{
#!   ignored
#! }
#! }
#! \details{
#! Function \code{ramsort} is a front-end to a couple of single-threaded sorting algorithms
#! that have been carefully implemented to be fast with and without \code{NA}s.
#! \cr
#! The default is a mergesort algorithm without copying (Sedgewick 8.4) for integer and double data
#! which requires 2x the RAM of its input vector (character or complex data are not supported).
#! Mergesort is fast, stable with a reliable runtime.
#! \cr
#! For integer data longer than a certain length we improve on mergesort by using a faster LSD
#! radixsort algorithm (Sedgewick 10.5) that uses 2x the RAM of its input vector plus 65536+1 integers.
#! \cr
#! For booleans, logicals, integers at or below the resolution of smallint and for factors below a certain number of levels
#! we use a key-index sort instead of mergesort or radix sort
#! (note that R has a (slower) key-index sort in \code{\link{sort.list}} available with confusingly named \code{method='radix'}
#! but the standard \code{\link{sort}} does not leverage it for factors (2-11.1).
#! If you call \code{keysort} directly, you should provide a known 'keyrange' directly to obtain the full speed.
#! \cr
#! Finally the user can request a sort method that minimizes memory use at the price of longer computation time
#! with \code{optimize='memory'} -- currently a shellsort.
#! }
#! \value{
#! integer scalar with the number of NAs. This is always 0 with has.na=FALSE
#! }
#! \references{
#! Robert Sedgewick (1997). Algorithms in C, Third edition. Addison-Wesley.
#! }
#! \author{
#! Jens Oehlschlägel
#! }
#! \note{
#! This function is called for its side-effects and breaks the functional programming paradigm. Use with care.
#! }
#!
#! \seealso{
#!   \code{\link{sort}}, \code{\link{ffsort}}, \code{\link{dfsort}}, \code{\link{ramorder}}
#! }
#! \examples{
#!    n <- 50
#!    x <- sample(c(NA, NA, 1:26), n, TRUE)
#!    sort(x)
#!    ramsort(x)
#!    x
#!
#!    \dontrun{
#!       message("Note how the datatype influences sorting speed")
#!       n <- 5e6
#!       x <- sample(1:26, n, TRUE)
#!
#!       y <- as.double(x)
#!       system.time(ramsort(y))
#!
#!       y <- as.integer(x)
#!       system.time(ramsort(y))
#!
#!       y <- as.short(x)
#!       system.time(ramsort(y))
#!
#!       y <- as.factor(letters)[x]
#!       system.time(ramsort(y))
#!    }
#! }
#! \keyword{univar}
#! \keyword{manip}
#! \keyword{arith}



#! \name{ramorder.default}
#! \alias{ramorder.default}
#! \alias{mergeorder.default}
#! \alias{radixorder.default}
#! \alias{keyorder.default}
#! \alias{shellorder.default}
#! \title{
#! Sorting: order R vector in-RAM and in-place
#! }
#! \description{
#! Function \code{ramorder} will order the input vector in-place (without making a copy) and return the number of NAs found
#! }
#! \usage{
#! \method{ramorder}{default}(x, i, has.na = TRUE, na.last = TRUE, decreasing = FALSE
#! , stable = TRUE, optimize = c("time", "memory"), VERBOSE = FALSE, \dots)
#! \method{mergeorder}{default}(x, i, has.na = TRUE, na.last = TRUE, decreasing = FALSE, \dots)
#! \method{radixorder}{default}(x, i, has.na = TRUE, na.last = TRUE, decreasing = FALSE, \dots)
#! \method{keyorder}{default}(x, i, keyrange=range(x, na.rm=has.na), has.na = TRUE, na.last = TRUE
#! , decreasing = FALSE, \dots)
#! \method{shellorder}{default}(x, i, has.na = TRUE, na.last = TRUE, decreasing = FALSE, stabilize=FALSE, \dots)
#! }
#! \arguments{
#!   \item{x}{
#! an atomic R vector
#! }
#!   \item{i}{
#!  a integer vector with a permuation of positions in x (you risk a crash if you violate this)
#! }
#!   \item{keyrange}{
#! an integer vector with two values giving the smallest and largest possible value in x, note that you should give this explicitely for best performance, relying on the default needs one pass over the data to determine the range
#! }
#!   \item{has.na}{
#! boolean scalar telling ramorder whether the vector might contain \code{NA}s.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{na.last}{
#! boolean scalar telling ramorder whether to order \code{NA}s last or first.
#! \emph{Note} that 'boolean' means that there is no third option \code{NA} as in \code{\link{order}}
#! }
#!   \item{decreasing}{
#! boolean scalar telling ramorder whether to order increasing or decreasing
#! }
#!   \item{stable}{
#!  set to false if stable ordering is not needed (may enlarge the set of ordering methods considered)
#! }
#!   \item{optimize}{
#! by default ramorder optimizes for 'time' which requires more RAM, set to 'memory'
#! to minimize RAM requirements and sacrifice speed
#! }
#!   \item{VERBOSE}{
#!   cat some info about chosen method
#! }
#!   \item{stabilize}{
#!   Set to \code{TRUE} for stabilizíng the result of shellorder (for equal keys the order values will be sorted, this only works if \code{i=1:n})
#! to minimize RAM requirements and sacrifice speed
#! }
#!   \item{\dots}{
#!   ignored
#! }
#! }
#! \details{
#! Function \code{ramorder} is a front-end to a couple of single-threaded ordering algorithms
#! that have been carefully implemented to be fast with and without \code{NA}s.
#! \cr
#! The default is a mergeorder algorithm without copying (Sedgewick 8.4) for integer and double data
#! which requires 2x the RAM of its input vector (character or complex data are not supported).
#! Mergeorder is fast, stable with a reliable runtime.
#! \cr
#! For integer data longer than a certain length we improve on mergeorder by using a faster LSD
#! radixorder algorithm (Sedgewick 10.5) that uses 2x the RAM of its input vector plus 65536+1 integers.
#! \cr
#! For booleans, logicals, integers at or below the resolution of smallint and for factors below a certain number of levels
#! we use a key-index order instead of mergeorder or radix order
#! (note that R has a (slower) key-index order in \code{\link{sort.list}} available with confusingly named \code{method='radix'}
#! but the standard \code{\link{order}} does not leverage it for factors (2-11.1).
#! If you call \code{keyorder} directly, you should provide a known 'keyrange' directly to obtain the full speed.
#! \cr
#! Finally the user can request a order method that minimizes memory use at the price of longer computation time
#! with \code{optimize='memory'} -- currently a shellorder.
#! }
#! \value{
#! integer scalar with the number of NAs. This is always 0 with has.na=FALSE
#! }
#! \references{
#! Robert Sedgewick (1997). Algorithms in C, Third edition. Addison-Wesley.
#! }
#! \author{
#! Jens Oehlschlägel
#! }
#! \note{
#! This function is called for its side-effects and breaks the functional programming paradigm. Use with care.
#! }
#!
#! \seealso{
#!   \code{\link{order}}, \code{\link{fforder}}, \code{\link{dforder}}, \code{\link{ramsort}}
#! }
#! \examples{
#!    n <- 50
#!    x <- sample(c(NA, NA, 1:26), n, TRUE)
#!    order(x)
#!    i <- 1:n
#!    ramorder(x, i)
#!    i
#!    x[i]
#!
#!    \dontrun{
#!       message("Note how the datatype influences sorting speed")
#!       n <- 1e7
#!       x <- sample(1:26, n, TRUE)
#!
#!       y <- as.double(x)
#!       i <- 1:n
#!       system.time(ramorder(y, i))
#!
#!       y <- as.integer(x)
#!       i <- 1:n
#!       system.time(ramorder(y, i))
#!
#!       y <- as.short(x)
#!       i <- 1:n
#!       system.time(ramorder(y, i))
#!
#!       y <- factor(letters)[x]
#!       i <- 1:n
#!       system.time(ramorder(y, i))
#!    }
#! }
#! \keyword{univar}
#! \keyword{manip}
#! \keyword{arith}



# orderes i inplace such that x[i] is sorted (x not affected)
mergeorder.default <- function(x, i, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...){
  force(x)
  force(i)
  n <- length(x)
  if (!is.integer(i))
    stop("i must be integer")
  if (length(i)!=n)
    stop("lengths of x and i don't match")
  #r <- range(i)
  #if (is.na(r[1]) || r[1]<1 || r[2]>n)
  #  stop("i must be C-code positions from 1:n")
  .Call("mergeorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}

# sorts x inplace
mergesort.default <- function(x, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...){
  force(x)
  .Call("mergesort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}



# orderes i inplace such that x[i] is sorted (x not affected)
shellorder.default <- function(x, i, has.na=TRUE, na.last=TRUE, decreasing=FALSE
, stabilize=FALSE, ...)  # note this requires that i was 1:length(x)
{
  force(x)
  force(i)
  n <- length(x)
  if (!is.integer(i))
    stop("i must be integer")
  if (length(i)!=n)
    stop("lengths of x and i don't match")
  #r <- range(i)
  #if (is.na(r[1]) || r[1]<1 || r[2]>n)
  #  stop("i must be C-code positions from 1:n")
  .Call("shellorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , stabilize  = as.logical(stabilize)
  , PACKAGE = "ff"
  )
}


# sorts x inplace
shellsort.default <- function(x, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...){
  force(x)
  .Call("shellsort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}

# orderes i inplace such that x[i] is sorted (x not affected)
# x must be integer with not too big range (<= 2^16 as a rule of thumb)
keyorder.default <- function(x, i, keyrange=range(x, na.rm=has.na)
, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...)
{
  force(x)
  force(i)
  n <- length(x)
  if (!is.integer(i))
    stop("i must be integer")
  if (length(i)!=n)
    stop("lengths of x and i don't match")
  #r <- range(i)
  #if (is.na(r[1]) || r[1]<1 || r[2]>n)
  #  stop("i must be C-code positions from 1:n")
  .Call("keyorder"
  , x = x
  , i = i
  , keyrange     = as.integer(keyrange)
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}

# sorts x inplace
# x must be integer with not too big range (<= 2^16 as a rule of thumb)
keysort.default <- function(x
, keyrange=range(x, na.rm=has.na)
, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...)
{
  force(x)
  .Call("keysort"
  , x = x
  , keyrange = as.integer(keyrange)
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}

# orderes i inplace such that x[i] is sorted (x not affected)
# x must be 4byte integer
radixorder.default <- function(x, i, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...){
  force(x)
  force(i)
  n <- length(x)
  if (!is.integer(i))
    stop("i must be integer")
  if (length(i)!=n)
    stop("lengths of x and i don't match")
  #r <- range(i)
  #if (is.na(r[1]) || r[1]<1 || r[2]>n)
  #  stop("i must be C-code positions from ^1:n")
  .Call("radixorder"
  , x = x
  , i = i
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}

# sorts x inplace
# x must be 4byte integer
radixsort.default <- function(x
, has.na=TRUE, na.last=TRUE, decreasing=FALSE, ...){
  force(x)
  .Call("radixsort"
  , x = x
  , has_na     = as.logical(has.na)
  , na_last    = as.logical(na.last)
  , decreasing = as.logical(decreasing)
  , PACKAGE = "ff"
  )
}


# sorts x inplace
# x must be atomic and currently without names
# returns the number of NAs encountered
ramsort.default <- function(x
, has.na=TRUE, na.last=TRUE, decreasing=FALSE
, optimize = c("time","memory")
, VERBOSE = FALSE
, ...
){
  n <- length(x)
  if (n){
    if (is.null(names(x))){
      optimize <- match.arg(optimize)
      if (optimize=="time"){
        v <- vmode(x)
        if (!is.na(.vNA[v]))
          has.na <- FALSE
        isint <- .ffmode[v]<=.ffmode["integer"]

        if (isint){
          l <- levels(x)
          if (is.null(l)){
            if (.ffmode[v]<=.ffmode["ushort"]){
              keyrange <- c(.vmin[v],.vmax[v])
              k <- diff(keyrange)+1L
            }else
              k <- 0L
          }else{
            k <- length(l)
            keyrange <- c(1L, k)
          }
          lk <- log(k,2)
          if (n>6300){
            # radixsort
            if ( k>0 && (lk<=17 || (lk<=26 && n > (3.742237e+05 + k*9.063064e-01)) )){
              if (VERBOSE)
                cat("ramsort selected keysort with", k ,"levels instead of radixsort\n")
              keysort(x, keyrange=keyrange
              , has.na=has.na, na.last=na.last, decreasing=decreasing)
            }else{
              if (VERBOSE)
                cat("ramsort selected radixsort\n")
              radixsort(x
              , has.na=has.na, na.last=na.last, decreasing=decreasing)
            }
          }else{
            # mergesort
            if (k>0 && lk<=26 && n> (-17.36857716 + k* 0.04508151) ){
              if (VERBOSE)
                cat("ramsort selected keysort with", k ,"levels instead of mergesort\n")
              keysort(x, keyrange=keyrange
              , has.na=has.na, na.last=na.last, decreasing=decreasing)
            }else{
              if (VERBOSE)
                cat("ramsort selected mergesort\n")
              mergesort(x
              , has.na=has.na, na.last=na.last, decreasing=decreasing)
            }
          }
        }else{
          if (VERBOSE)
            cat("ramorder selected mergesort\n")
          mergesort(x
          , has.na=has.na, na.last=na.last, decreasing=decreasing)
        }
      }else{
        shellsort(x
        , has.na=has.na, na.last=na.last, decreasing=decreasing)
      }
    }else{
      stop("handling names not implemented")
    }
  }
}


# orderes i inplace such that x[i] is sorted (x not affected)
# x must be atomic and currently without names
# i must be integer with valid positions (R counting 1..n)
# returns the number of NAs encountered
ramorder.default <- function(x, i
, has.na=TRUE, na.last=TRUE, decreasing=FALSE
, stable = TRUE
, optimize = c("time","memory")
, VERBOSE = FALSE
, ...
){
  n <- length(x)
  if (n){
    optimize <- match.arg(optimize)
    if (stable || optimize=="time"){
      v <- vmode(x)
      if (!is.na(.vNA[v]))
        has.na <- FALSE
      isint <- .ffmode[v]<=.ffmode["integer"]
      if (isint){
        l <- levels(x)
        if (is.null(l)){
          if (.ffmode[v]<=.ffmode["ushort"]){
            keyrange <- c(.vmin[v],.vmax[v])
            k <- diff(keyrange)+1L
          }else
            k <- 0L
        }else{
          k <- length(l)
          keyrange <- c(1L, k)
        }
        lk <- log(k,2)
        if (n>6000){
          # radixorder
          if ( k>0 && (lk<=17 || (lk<=26 && n > 2^(-250.32441 - 10.80749*lk + 108.44324*sqrt(lk)) ) ) ){
            if (VERBOSE)
              cat("ramorder selected keyorder with", k ,"levels instead of radixorder\n")
            keyorder(x, i, keyrange=keyrange
            , has.na=has.na, na.last=na.last, decreasing=decreasing)
          }else{
            if (VERBOSE)
              cat("ramorder selected radixorder\n")
            radixorder(x, i
            , has.na=has.na, na.last=na.last, decreasing=decreasing)
          }
        }else{
          # mergeorder
          if (k>0 && lk<=26 && n> 2^(7.17728849 - 0.43158315*lk + 0.04323707*lk^2) ){
            if (VERBOSE)
              cat("ramorder selected keyorder with", k ,"levels instead of mergeorder\n")
            keyorder(x, i, keyrange=keyrange
            , has.na=has.na, na.last=na.last, decreasing=decreasing)
          }else{
            if (VERBOSE)
              cat("ramorder selected mergeorder\n")
            mergeorder(x, i
            , has.na=has.na, na.last=na.last, decreasing=decreasing)
          }
        }
      }else{
        if (VERBOSE)
          cat("ramorder selected mergeorder\n")
        mergeorder(x, i
        , has.na=has.na, na.last=na.last, decreasing=decreasing)
      }
    }else{
      if (VERBOSE)
        cat("ramorder selected shellorder\n")
      shellorder(x, i
      , has.na=has.na, na.last=na.last, decreasing=decreasing)
    }
  }
}




#! \name{ffsort}
#! \alias{ffsort}
#! \title{
#! Sorting of ff vectors
#! }
#! \description{
#! Sorting: sort an ff vector -- optionally in-place
#! }
#! \usage{
#! ffsort(x
#! , aux = NULL
#! , has.na = TRUE
#! , na.last = TRUE
#! , decreasing = FALSE
#! , inplace = FALSE
#! , decorate = FALSE
#! , BATCHBYTES = getOption("ffmaxbytes")
#! , VERBOSE = FALSE
#! )
#! }
#! \arguments{
#!   \item{x}{
#!     an ff vector
#! }
#!   \item{aux}{
#!     NULL or an ff vector of the same type for temporary storage
#! }
#!   \item{has.na}{
#! boolean scalar telling ffsort whether the vector might contain \code{NA}s.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{na.last}{
#! boolean scalar telling ffsort whether to sort \code{NA}s last or first.
#! \emph{Note} that 'boolean' means that there is no third option \code{NA} as in \code{\link{sort}}
#! }
#!   \item{decreasing}{
#! boolean scalar telling ffsort whether to sort increasing or decreasing
#! }
#!   \item{inplace}{
#! boolean scalar telling ffsort whether to sort the original ff vector (\code{TRUE})
#! or to create a sorted copy (\code{FALSE}, the default)
#! }
#!   \item{decorate}{
#! boolean scalar telling ffsort whether to decorate the returned ff vector with \code{\link{is.sorted}}
#! and \code{\link{na.count}} attributes.
#! }
#!   \item{BATCHBYTES}{
#!   maximum number of RAM bytes ffsort should try not to exceed
#! }
#!   \item{VERBOSE}{
#!   cat some info about the sorting
#! }
#! }
#! \details{
#!   ffsort tries to sort the vector in-RAM respecting the BATCHBYTES limit.
#!   If a fast sort it not possible, it uses a slower in-place sort (shellsort).
#!   If in-RAM is not possible, it uses (a yet simple) out-of-memory algorithm.
#!   Like \code{\link{ramsort}} the in-RAM sorting method is choosen depending on context information.
#!   If a key-index sort can be used, ffsort completely avoids merging disk based subsorts.
#!   If argument \code{decorate=TRUE} is used, then \code{na.count(x)} will return the number of NAs
#!   and \code{is.sorted(x)} will return TRUE if the sort was done with \code{na.last=TRUE} and \code{decreasing=FALSE}.
#! }
#! \note{
#!   the ff vector may not have a names attribute
#! }
#! \value{
#!   An ff vector -- optionally decorated with \code{\link{is.sorted}} and \code{\link{na.count}}, see argument 'decorate'
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#!
#! \seealso{
#!   \code{\link{ramsort}}, \code{\link{fforder}}, \code{\link{ffdfsort}}
#! }
#! \examples{
#!    n <- 1e6
#!    x <- ff(c(NA, 999999:1), vmode="double", length=n)
#!    x <- ffsort(x)
#!    x
#!    is.sorted(x)
#!    na.count(x)
#!    x <- ffsort(x, decorate=TRUE)
#!    is.sorted(x)
#!    na.count(x)
#!    x <- ffsort(x, BATCHBYTES=n, VERBOSE=TRUE)
#! }
#! \keyword{univar}
#! \keyword{manip}
#! \keyword{arith}
#! \keyword{ IO }
#! \keyword{ data }

# returns the sorted ff vector optionally with decorations 'na.count' and 'is.sorted' (the latter only if na.last=TRUE and decreasing=FALSE)
ffsort <- function(
  x
, aux = NULL
, has.na=TRUE, na.last=TRUE, decreasing=FALSE
, inplace = FALSE
, decorate = FALSE
, BATCHBYTES  = getOption("ffmaxbytes")
, VERBOSE     = FALSE
){
  stopifnot(is.null(vw(x)))
  stopifnot(is.null(dim(x)))
  n <- length(x)
  v <- vmode(x)
  if (!is.na(.vNA[v]))
    has.na <- FALSE

  isint <- .ffmode[v]<=.ffmode["integer"]

  if (isint){
    l <- levels(x)
    if (is.null(l)){
      if (.ffmode[v]<=.ffmode["ushort"]){
        keyrange <- c(.vmin[v],.vmax[v])
        k <- diff(keyrange)+1L
      }else
        k <- 0L
    }else{
      k <- length(l)
      if (is.na(.vNA[v]))
        keyrange <- c(1L, k)
      else
        keyrange <- c(0L, k-1L)
    }
  }else{
    k <- 0L
  }

  if (is.null(names(x))){

    recbytes <- .rambytes[v]

    if (k && (k < .vvalues["short"])){
	  method <- 3L
      keybytes <- (k+2L) * .rambytes["integer"]
      maxbytes <- max(BATCHBYTES, 3L*keybytes)
      maxordersize <- (maxbytes - 2L*keybytes) %/% recbytes
      if (maxordersize>n)
        BATCHSIZE <- bbatch(n, maxordersize)$b
      else
        BATCHSIZE <- n
    }else{
      keybytes <- (.vvalues["short"]+1L) * .rambytes["integer"]
      maxbytes <- max(BATCHBYTES, 3L*keybytes)
      allbytes <- 2 * n * recbytes
      if (isint && allbytes + keybytes <= maxbytes){
        method <- 2L  # radixsort
        BATCHSIZE <- n
      }else if(allbytes <= maxbytes){
        method <- 0L  # mergesort
        BATCHSIZE <- n
      }else{
        method <- 1L  # shellsort
        BATCHSIZE <- ceiling(n / 2 ^ ceiling(max(0, log(n*(recbytes/maxbytes), 2))))
      }
    }
    if (BATCHSIZE>=n){
     # tuning: do it in-ram
     y <- read.ff(x, 1L, n)
      if (VERBOSE)
        cat("method=")
      nNA <- ramsort(y, has.na = has.na, na.last = na.last, decreasing = decreasing, optimize = if (method==1L) "memory" else "time", VERBOSE = VERBOSE)
      if (!inplace)
        x <- clone(x, initdata=NULL)
      write.ff(x, 1L, y)
      rm(y)
    }else{
      # fallback: do it on-disk
      if (!inplace)
        x <- clone(x)
			else
				open(x, assert=TRUE)

        if (k && (k < .vvalues["short"])){

          if (VERBOSE)
            cat("method=ffkeysort  BATCHSIZE=", BATCHSIZE, "\n", sep="")
          nNA <- .Call("ffkeysort"
          , .ffmode[v]
          , ff_           = attr(x, "physical")
          , left_         = 1L
          , right_        = n
          , keyrange_     = as.integer(keyrange)
          , ordersize_    = as.integer(BATCHSIZE)
          , has_na_       = as.logical(has.na)
          , na_last_      = as.logical(na.last)
          , decreasing_   = as.logical(decreasing)
          , PACKAGE="ff"
          )

        }else{

          if (VERBOSE)
            cat("method=", c("merge","shell","radix","key","quick")[method+1L],"  BATCHSIZE=", BATCHSIZE, "\n", sep="")
          if (is.null(aux)){
            aux <- clone(x, pattern="tmpffsort")
          }else{
            if (vmode(aux) != v)
            stop("vmode(aux) does not match vmode(x)")
            if (length(aux) != n)
            stop("length(aux) does not match length(x)")
						open(aux, assert=TRUE)
          }
          nNA <- .Call("ffsortmerge"
          , ffmode_       = .ffmode[v]
          , ff_           = attr(x, "physical")
          , auxff_        = attr(aux, "physical")
          , left_         = 1L
          , right_        = n
          , method_       = as.integer(method)
          , keyrange_     = NULL  # in sorting, we always use ffkeysort if we can use keysort
          , ordersize_    = as.integer(BATCHSIZE)  # here we use as much as possible
          , mergesize_    = pagesize(x)            # while here we save as much as possible to support the file system cache
          , has_na_       = as.logical(has.na)
          , na_last_      = as.logical(na.last)
          , decreasing_   = as.logical(decreasing)
          , PACKAGE="ff"
          )

        }
    }

  }else{
    stop("handling names not implemented")
  }

  if (decorate){
    na.count(x) <- nNA
    is.sorted(x) <- na.last && !decreasing
  }
  x
}


#! \name{fforder}
#! \alias{fforder}
#! \title{
#! Sorting: order from ff vectors
#! }
#! \description{
#! Returns order with regard to one or more ff vectors
#! }
#! \usage{
#! fforder(...
#! , index = NULL
#! , use.index = NULL
#! , aux = NULL
#! , auxindex = NULL
#! , has.na = TRUE
#! , na.last = TRUE
#! , decreasing = FALSE
#! , BATCHBYTES = getOption("ffmaxbytes")
#! , VERBOSE = FALSE
#! )
#! }
#! \arguments{
#!   \item{\dots}{
#!   one of more ff vectors which define the order
#! }
#!   \item{index}{
#!   an optional ff integer vector used to store the order output
#! }
#!   \item{use.index}{
#!   A boolean flag telling fforder whether to use the positions in 'index' as input.
#!   If you do this, it is your responsibility to assure legal positions - otherwise you risk a crash.
#! }
#!   \item{aux}{
#!   An optional named list of ff vectors that can be used for temporary copying
#!   -- the names of the list identify the \code{\link{vmode}s} for which the respective ff vector is suitable.
#! }
#!   \item{auxindex}{
#!   An optional ff intger vector for temporary storage of integer positions.
#! }
#!   \item{has.na}{
#! boolean scalar telling fforder whether the vector might contain \code{NA}s.
#! \emph{Note} that you risk a crash if there are unexpected \code{NA}s with \code{has.na=FALSE}
#! }
#!   \item{na.last}{
#! boolean scalar telling fforder whether to order \code{NA}s last or first.
#! \emph{Note} that 'boolean' means that there is no third option \code{NA} as in \code{\link{order}}
#! }
#!   \item{decreasing}{
#! boolean scalar telling fforder whether to order increasing or decreasing
#! }
#!   \item{BATCHBYTES}{
#!   maximum number of RAM bytes fforder should try not to exceed
#! }
#!   \item{VERBOSE}{
#!   cat some info about the ordering
#! }
#! }
#! \details{
#!   fforder tries to order the vector in-RAM, if not possible it uses (a yet simple) out-of-memory algorithm.
#!   Like \code{\link{ramorder}} the in-RAM ordering method is choosen depending on context information.
#! }
#! \value{
#!   An ff vector with the positions that ore required to sort the input as specified
#! -- with an attribute \code{\link{na.count}} with as many values as columns in \dots
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#!
#! \seealso{
#!   \code{\link{ramorder}}, \code{\link{ffsort}}, \code{\link{ffdforder}}, \code{\link{ffindexget}}
#! }
#! \examples{
#!    x <- ff(sample(1e5, 1e6, TRUE))
#!    y <- ff(sample(1e5, 1e6, TRUE))
#!    d <- ffdf(x, y)
#!
#!    i <- fforder(y)
#!    y[i]
#!    i <- fforder(x, index=i)
#!    x[i]
#!    d[i,]
#!
#!    i <- fforder(x, y)
#!    d[i,]
#!
#!    i <- ffdforder(d)
#!    d[i,]
#!
#!    rm(x, y, d, i)
#!    gc()
#! }
#! \keyword{univar}
#! \keyword{manip}
#! \keyword{arith}
#! \keyword{ IO }
#! \keyword{ data }


# the returned vector carries an attribute 'na.count' with as many values as columns in ...
fforder <- function(
  ...
, index = NULL
, use.index = NULL
, aux = NULL
, auxindex = NULL
, has.na=TRUE, na.last=TRUE, decreasing=FALSE
, BATCHBYTES  = getOption("ffmaxbytes")
, VERBOSE     = FALSE
)
{
  if (is.null(use.index))
    use.index <- !is.null(index)

  l <- list(...)
  nl <- length(l)
  if (is.null(names(l)))
    names(l) <- paste("<", 1:nl, ">", sep="")
  if (nl==1){
    if (!is.null(aux) && !is.list(aux)){
      aux <- list(aux)
      names(aux) <- vmode(l[[1]])
    }
  }else{
    has.na <- rep(has.na, length.out=nl)
    na.last <- rep(na.last, length.out=nl)
    decreasing <- rep(decreasing, length.out=nl)
  }

  analyze <- function(x, has.na){
    stopifnot(is.null(vw(x)))
    stopifnot(is.null(dim(x)))
    n <- length(x)
    v <- vmode(x)
    isint <- .ffmode[v]<=.ffmode["integer"]

    if (isint){
      l <- levels(x)
      if (is.null(l)){
        if (.ffmode[v]<=.ffmode["ushort"]){
          keyrange <- c(.vmin[v],.vmax[v])
          k <- diff(keyrange)+1L
        }else{
          k <- 0L
          keyrange <- NULL
        }
      }else{
        k <- length(l)
        if (is.na(.vNA[v]))
          keyrange <- c(1L, k)
        else
          keyrange <- c(0L, k-1L)
      }
    }else{
      k <- 0L
      keyrange <- NULL
    }

    recvalbytes <- .rambytes[v]
    recindbytes <- .rambytes["integer"]

    if (k && (k < .vvalues["short"])){
      method <- 3L # keyorder
      keybytes <- (k+2L) * .rambytes["integer"]
      maxbytes <- max(BATCHBYTES, 3L*keybytes)
      allbytes <- n*(recvalbytes+2*recindbytes) + keybytes
    }else if (isint){
      method <- 2L  # radixorder
      keybytes <- (.vvalues["short"]+1L) * .rambytes["integer"]
      maxbytes <- max(BATCHBYTES, 3L*keybytes)
      allbytes <- n*(recvalbytes+2*recindbytes) + keybytes
    }else{
      method <- 0L  # mergeorder
      maxbytes <- BATCHBYTES
      allbytes <- n*(recvalbytes+2*recindbytes)
    }

    pieces <- 2 ^ ceiling(max(0, log(allbytes/maxbytes, 2)))
    if (pieces>1 && !use.index){
      maxbytes <- BATCHBYTES
      allbytes <- as.double(n)*(recvalbytes+recindbytes)
      shellpieces <- 2 ^ ceiling(max(0, log(allbytes/maxbytes, 2)))
      if (shellpieces<pieces){
        method <- 1L
        pieces <- shellpieces
      }
    }
    BATCHSIZE <- ceiling(n / pieces)

    list(v=v, n=n, method=method, keyrange=keyrange, BATCHSIZE=BATCHSIZE)
  }

  a <- do.call("cbind", lapply(l, analyze))
  n <- unique(unlist(a["n",]))
  stopifnot(length(n)==1)

  if (is.null(index)){
    index <- clone(l[[1]], vmode="integer", initdata=NULL)
	levels(index) <- NULL
    use.index <- FALSE
  }else{
    stopifnot(vmode(index)=="integer")
    stopifnot(length(index)==n)
  }

  maxBATCHSIZE <- max(unlist(a["BATCHSIZE", ]))
  nNA <- integer(nl)

  if (maxBATCHSIZE==n){
    # tuning: do it in-ram

    for (i in nl:1){
      x <- read.ff(l[[i]], 1L, n)
      if (i==nl){
        if (use.index){
          o <- read.ff(index, 1L, n)
        }else{
          o <- 1:n
        }
      }

      if (VERBOSE){
        BATCHSIZE <- a[["BATCHSIZE", i]]
        cat("column=", names(l)[i], "  method=ramorder  BATCHSIZE=", BATCHSIZE, "\n", sep="")
      }
      nNA[i] <- ramorder(x, o
      , has.na = has.na[i]
      , na.last = na.last[i]
      , decreasing = decreasing[i]
      , stable = TRUE
      , optimize = "time"
      )
    }
    write.ff(index, 1L, o)

  }else{
    # fallback: do it on-disk

    if (is.null(auxindex)){
      auxindex <- clone(index, pattern="tmpfforder")
    }else{
      if (vmode(auxindex) != "integer")
        stop("vmode(auxindex) does not match 'integer'")
      if (length(auxindex) != n)
        stop("length(auxindex) does not match length(x)")
			open(index, assert=TRUE)
			open(auxindex, assert=TRUE)
    }

    for (i in nl:1){
      x <- clone(l[[i]])
      if (i<nl){
        x <- ffindexget(x, index, BATCHSIZE=BATCHSIZE, VERBOSE=VERBOSE)
        use.index <- TRUE
      }
      BATCHSIZE <- a[["BATCHSIZE", i]]
      method <- a[["method", i]]
      v <- a[["v", i]]
      if (!is.na(.vNA[v]))
        has.na[i] <- FALSE

      if (BATCHSIZE<n){
        if (is.null(aux[[v]])){
          aux[[v]] <- clone(x, pattern="tmpfforder")
        }else{
          if (length(aux[[v]]) != n)
            stop("length(aux) does not match length(x)")
					open(aux[[v]], assert=TRUE)
        }
      }

      if (VERBOSE)
        cat("column=", names(l)[i], "  method=", c("merge","shell","radix","key","quick")[method+1L],"  BATCHSIZE=", BATCHSIZE, "\n", sep="")
      nNA[i] <- .Call("ffordermerge"
      , ffmode_       = .ffmode[v]
      , ff_           = attr(x, "physical")
      , index_        = attr(index, "physical")
      , auxff_        = attr(aux[[v]], "physical")
      , auxindex_     = attr(auxindex, "physical")
      , left_         = 1L
      , right_        = n
      , method_       = as.integer(method)
      , keyrange_     = a[["keyrange", i]]
      , ordersize_    = as.integer(BATCHSIZE)  # here we use as much as possible
      , mergesize_    = pagesize(x)            # while here we save as much as possible to support the file system cache
      , orderindex_   = as.logical(use.index)
      , has_na_       = as.logical(has.na[i])
      , na_last_      = as.logical(na.last[i])
      , decreasing_   = as.logical(decreasing[i])
      , PACKAGE="ff"
      )
    }

  }
  attr(index, "na.count") <- nNA
  index
}





#! \name{ffindexorder}
#! \alias{ffindexorder}
#! \alias{ffindexordersize}
#! \title{
#!   Sorting: chunked ordering of integer suscript positions
#! }
#! \description{
#!   Function \code{ffindexorder} will calculate chunkwise the order positions to sort all positions in a chunk ascending.
#!   \cr
#!   Function \code{ffindexordersize} does the calculation of the chunksize for \code{ffindexorder}.
#! }
#! \usage{
#! ffindexordersize(length, vmode, BATCHBYTES = getOption("ffmaxbytes"))
#! ffindexorder(index, BATCHSIZE, FF_RETURN = NULL, VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{index}{
#!   A \code{\link{ff}} integer vector with integer subscripts.
#! }
#!   \item{BATCHSIZE}{
#!   Limit for the chunksize (see details)
#! }
#!   \item{BATCHBYTES}{
#!   Limit for the number of bytes per batch
#! }
#!   \item{FF_RETURN}{
#!   Optionally an \code{\link{ff}} integer vector in which the chunkwise order positions are stored.
#! }
#!   \item{VERBOSE}{
#!   Logical scalar for activating verbosing.
#! }
#!   \item{length}{
#!   Number of elements in the index
#! }
#!   \item{vmode}{
#!   The \code{\link{vmode}} of the ff vector to which the index shall be applied with \code{\link{ffindexget}} or \code{\link{ffindexset}}
#! }
#! }
#! \details{
#!   Accessing integer positions in an ff vector is a non-trivial task, because it could easily lead to random-access to a disk file.
#!   We avoid random access by loading batches of the subscript values into RAM, order them ascending, and only then access the ff values on disk.
#!   Such an ordering can be done on-the-fly by \code{\link{ffindexget}} or it can be created upfront with \code{ffindexorder}, stored and re-used,
#!   similar to storing and using hybrid index information with \code{\link{as.hi}}.
#! }
#! \value{
#!   Function \code{ffindexorder} returns an ff integer vector with an attribute \code{BATCHSIZE} (the chunksize finally used, not the one given with argument \code{BATCHSIZE}).
#!   \cr
#!   Function \code{ffindexordersize} returns a balanced batchsize as returned from \code{\link{bbatch}}.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!    \code{\link{ffindexget}},  \code{\link{as.hi}},  \code{\link{bbatch}}
#! }
#! \examples{
#!      x <- ff(sample(40))
#!      message("fforder requires sorting")
#!      i <- fforder(x)
#!      message("applying this order i is done by ffindexget")
#!      x[i]
#!      message("applying this order i requires random access, 
#!        therefore ffindexget does chunkwise sorting")
#!      ffindexget(x, i)
#!      message("if we want to apply the order i multiple times,
#!        we can do the chunkwise sorting once and store it")
#!      s <- ffindexordersize(length(i), vmode(i), BATCHBYTES = 100)
#!      o <- ffindexorder(i, s$b)
#!      message("this is how the stored chunkwise sorting is used")
#!      ffindexget(x, i, o)
#!      message("")
#!      rm(x,i,s,o)
#!      gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


ffindexordersize <- function(length, vmode, BATCHBYTES = getOption("ffmaxbytes")){
  recvalbytes <- .rambytes[vmode]
  recindbytes <- .rambytes["integer"]
  bbatch(length, BATCHBYTES / (recvalbytes+2*recindbytes))  #assuming chunkorder method 'quick'
}

ffindexorder <- function(
  index
, BATCHSIZE
, FF_RETURN = NULL
, VERBOSE     = FALSE
){
  if (VERBOSE)
    starttime <- proc.time()[3]

  stopifnot(is.null(vw(index)))
  stopifnot(vmode(index)=="integer")
  n <- length(index)


	open(index, assert=TRUE)
  if (is.null(FF_RETURN)){
    FF_RETURN <- clone(index, vmode="integer", initdata=NULL)
  }else{
    stopifnot(is.null(vw(FF_RETURN)))
    stopifnot(vmode(FF_RETURN)=="integer")
    stopifnot(length(FF_RETURN)==n)
		open(FF_RETURN, assert=TRUE)
  }

  .Call("ffchunkorder"
  , index_        = attr(index, "physical")
  , indexorder_   = attr(FF_RETURN, "physical")
  , indexsize_    = n
  , method_       = 4L  # quicksort
  , ordersize_    = as.integer(BATCHSIZE)
  , PACKAGE="ff"
  )
  attr(FF_RETURN, "BATCHSIZE") <- as.integer(BATCHSIZE)

  if (VERBOSE)
    cat("ffindexorder", proc.time()[3] - starttime, "seconds\n")

  FF_RETURN
}


#! \name{ffindexget}
#! \alias{ffindexget}
#! \alias{ffindexset}
#! \title{
#!   Reading and writing ff vectors using ff subscripts
#! }
#! \description{
#!   Function \code{ffindexget} allows to extract elements from an ff vector according to positive integer suscripts stored in an ff vector.
#!   \cr
#!   Function \code{ffindexset} allows the inverse operation: assigning to elements of an ff vector according to positive integer suscripts stored in an ff vector.
#!   These functions allow more control than the method dispatch of \code{[} and  \code{[<-} if an ff integer subscript is used.
#! }
#! \usage{
#! ffindexget(x, index, indexorder = NULL, FF_RETURN = NULL
#! , BATCHSIZE = NULL, BATCHBYTES = getOption("ffmaxbytes"), VERBOSE = FALSE)
#! ffindexset(x, index, value, indexorder = NULL
#! , BATCHSIZE = NULL, BATCHBYTES = getOption("ffmaxbytes"), VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{x}{
#!   A \code{\link{ff}} vector containing the elements
#! }
#!   \item{index}{
#!   A \code{\link{ff}} integer vector with integer subscripts in the range from \code{1} to \code{length(x)}.
#! }
#!   \item{value}{
#!   An \code{\link{ff}} vector of the same \code{\link{vmode}} as x containing the values to be assigned
#! }
#!   \item{indexorder}{
#!   Optionally the return value of \code{\link{ffindexorder}}, see details
#! }
#!   \item{FF_RETURN}{
#!   Optionally an \code{\link{ff}} vector of the same \code{\link{vmode}} as x in which the returned values shall be stored, see details.
#! }
#!   \item{BATCHSIZE}{
#!   Optinal limit for the batchsize (see details)
#! }
#!   \item{BATCHBYTES}{
#!   Limit for the number of bytes per batch
#! }
#!   \item{VERBOSE}{
#!   Logical scalar for verbosing
#! }
#! }
#! \details{
#!   Accessing integer positions in an ff vector is a non-trivial task, because it could easily lead to random-access to a disk file.
#!   We avoid random access by loading batches of the subscript values into RAM, order them ascending, and only then access the ff values on disk.
#!   Since ordering is expensive, it may pay to do the batched ordering once upfront and then re-use it with  \code{\link{ffindexorder}},
#!   similar to storing and using hybrid index information with \code{\link{as.hi}}.
#! }
#! \value{
#!   Function \code{ffindexget} returns an ff vector with the extracted elements.
#!   \cr
#!   Function \code{ffindexset} returns the ff vector in which we have updated values.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{Extract.ff}}, \code{\link{ffdfindexget}}, \code{\link{ffindexorder}}
#! }
#! \examples{
#! message("ff integer subscripts with ff return/assign values")
#! x <- ff(factor(letters))
#! i <- ff(2:9)
#! xi <- x[i]
#! xi
#! xi[] <- NA
#! xi
#! x[i] <- xi
#! x
#! message("ff integer subscripts: more control with ffindexget/ffindexset")
#! xi <- ffindexget(x, i, FF_RETURN=xi)
#! x <- ffindexset(x, i, xi)
#! rm(x, i, xi)
#! gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

ffindexget <- function(
  x
, index
, indexorder = NULL
, FF_RETURN = NULL
, BATCHSIZE   = NULL
, BATCHBYTES  = getOption("ffmaxbytes")
, VERBOSE     = FALSE
){

  stopifnot(is.null(vw(x)))
  stopifnot(is.null(dim(x)))
  stopifnot(is.null(names(x)))

  n <- length(index)
  N <- length(x)
  stopifnot(vmode(index)=="integer")
  v <- vmode(x)

  if (is.null(FF_RETURN)){
    FF_RETURN <- clone(x, length=n, initdata=NULL)
  }else{
    stopifnot(vmode(FF_RETURN)==v)
		open(FF_RETURN, assert=TRUE)
    if (length(FF_RETURN)!=n)
      length(FF_RETURN) <- n
  }

  recvalbytes <- .rambytes[v]
  recindbytes <- .rambytes["integer"]
  theobytes.inram <- 2*recvalbytes+2*recindbytes

  if ( n*theobytes.inram<=BATCHBYTES && n>=N && (is.null(BATCHSIZE) || BATCHSIZE>=n) ){
    # tuning: do it in-ram
    if (VERBOSE)
      cat("method=inram BATCHSIZE=", n, "\n", sep="")

    #initram <- memory.size()

    i <- read.ff(index, 1L, n)
    #cat("i "); print(memory.size()-initram)

    r <- read.ff(x, 1L, N)
    #cat("r+i "); print(memory.size()-initram)

    r <- r[i]
    #cat("r[i] "); print(memory.size()-initram)

    rm(i)

    FF_RETURN <- write.ff(FF_RETURN, 1L, r, add = FALSE)
    #cat("after write "); print(memory.size()-initram)

  }else{

    # fallback: do it on-disk
    if (is.null(indexorder)){
      BATCHSIZE <- min(BATCHSIZE, ffindexordersize(length=n, vmode=v, BATCHBYTES=BATCHBYTES)$b)
      if (VERBOSE)
        cat("method=", c("merge","shell","radix","key","quick")[5L],"  BATCHSIZE=", BATCHSIZE, "\n", sep="")
    }else{
      stopifnot(is.null(vw(indexorder)))
      stopifnot(vmode(indexorder)=="integer")
      stopifnot(length(indexorder)==n)
      BATCHSIZE <- attr(indexorder, "BATCHSIZE")
      if (VERBOSE)
        cat("method=indexorder  BATCHSIZE=", BATCHSIZE, "\n", sep="")
				open(indexorder, assert=TRUE)
    }

		open(x, assert=TRUE)
		open(index, assert=TRUE)
    .Call("ffindexget"
    , ffmode_       = .ffmode[v]
    , baseff_       = attr(x, "physical")
    , returnff_     = attr(FF_RETURN, "physical")
    , index_        = attr(index, "physical")
    , auxindex_     = attr(indexorder, "physical")  # is NULL for NULL
    , offset        = 1L
    , left_         = 1L
    , right_        = n
    , method_       = 4L  # quicksort
    , ordersize_    = as.integer(BATCHSIZE)
    , PACKAGE="ff"
    )
  }

  FF_RETURN
}


ffindexset <- function(
  x
, index
, value
, indexorder = NULL
, BATCHSIZE   = NULL
, BATCHBYTES  = getOption("ffmaxbytes")
, VERBOSE     = FALSE
){

  stopifnot(is.null(vw(x)))
  stopifnot(is.null(dim(x)))
  n <- length(index)
  N <- length(x)
  stopifnot(vmode(index)=="integer")
  v <- vmode(x)

  stopifnot(is.null(vw(value)))
  stopifnot(vmode(value)==v)
  stopifnot(length(value)==n)

  recvalbytes <- .rambytes[v]
  recindbytes <- .rambytes["integer"]
  theobytes.inram <- 2*recvalbytes+2*recindbytes

  if ( n*theobytes.inram<=BATCHBYTES && n>=N && (is.null(BATCHSIZE) || BATCHSIZE>=n) ){
    # tuning: do it in-ram
    if (VERBOSE)
      cat("method=inram BATCHSIZE=", n, "\n", sep="")

    #initram <- memory.size()

    v <- read.ff(x, 1L, N)
    #cat("v "); print(memory.size()-initram)
    v[read.ff(index, 1L, n)] <- read.ff(value, 1L, n)
    #cat("v[index] <- value "); print(memory.size()-initram)

    write.ff(x, 1L, v)
    #cat("after write "); print(memory.size()-initram)

  }else{

    # fallback: do it on-disk
    if (is.null(indexorder)){
      BATCHSIZE <- min(BATCHSIZE, ffindexordersize(length=n, vmode=v, BATCHBYTES=BATCHBYTES)$b)
      if (VERBOSE)
        cat("method=", c("merge","shell","radix","key","quick")[5L],"  BATCHSIZE=", BATCHSIZE, "\n", sep="")
    }else{
      stopifnot(is.null(vw(indexorder)))
      stopifnot(vmode(indexorder)=="integer")
      stopifnot(length(indexorder)==n)
      BATCHSIZE <- attr(indexorder, "BATCHSIZE")
      if (VERBOSE)
        cat("method=indexorder  BATCHSIZE=", BATCHSIZE, "\n", sep="")
			open(indexorder, assert=TRUE)
    }
		open(x, assert=TRUE)
		open(value, assert=TRUE)
		open(index, assert=TRUE)

    .Call("ffindexset"
    , ffmode_       = .ffmode[v]
    , baseff_       = attr(x, "physical")
    , valueff_      = attr(value, "physical")
    , index_        = attr(index, "physical")
    , auxindex_     = attr(indexorder, "physical")  # is NULL for NULL
    , offset        = 1L
    , left_         = 1L
    , right_        = n
    , method_       = 4L  # quicksort
    , ordersize_    = as.integer(BATCHSIZE)
    , PACKAGE="ff"
    )
  }

  x
}



#! \name{ffdfindexget}
#! \alias{ffdfindexget}
#! \alias{ffdfindexset}
#! \title{
#!   Reading and writing ffdf data.frame using ff subscripts
#! }
#! \description{
#!   Function \code{ffdfindexget} allows to extract rows from an ffdf data.frame according to positive integer suscripts stored in an ff vector.
#!   \cr
#!   Function \code{ffdfindexset} allows the inverse operation: assigning to rows of an ffdf data.frame according to positive integer suscripts stored in an ff vector.
#!   These functions allow more control than the method dispatch of \code{[} and  \code{[<-} if an ff integer subscript is used.
#! }
#! \usage{
#!   ffdfindexget(x, index, indexorder = NULL, autoindexorder = 3, FF_RETURN = NULL
#!   , BATCHSIZE = NULL, BATCHBYTES = getOption("ffmaxbytes"), VERBOSE = FALSE)
#!   ffdfindexset(x, index, value, indexorder = NULL, autoindexorder = 3
#!   , BATCHSIZE = NULL, BATCHBYTES = getOption("ffmaxbytes"), VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{x}{
#!   A \code{\link{ffdf}} data.frame containing the elements
#! }
#!   \item{index}{
#!   A \code{\link{ff}} integer vector with integer subscripts in the range from \code{1} to \code{length(x)}.
#! }
#!   \item{value}{
#!   A \code{\link{ffdf}} data.frame like \code{x} with the rows to be assigned
#! }
#!   \item{indexorder}{
#!   Optionally the return value of \code{\link{ffindexorder}}, see details
#! }
#!   \item{autoindexorder}{
#!   The minimum number of columns (which need chunked indexordering) for which we switch from on-the-fly ordering to stored \code{\link{ffindexorder}}
#! }
#!   \item{FF_RETURN}{
#!   Optionally an \code{\link{ffdf}} data.frame of the same type as x in which the returned values shall be stored, see details.
#! }
#!   \item{BATCHSIZE}{
#!   Optinal limit for the batchsize (see details)
#! }
#!   \item{BATCHBYTES}{
#!   Limit for the number of bytes per batch
#! }
#!   \item{VERBOSE}{
#!   Logical scalar for verbosing
#! }
#! }
#! \details{
#!   Accessing rows of an ffdf data.frame identified by integer positions in an ff vector is a non-trivial task, because it could easily lead to random-access to disk files.
#!   We avoid random access by loading batches of the subscript values into RAM, order them ascending, and only then access the ff values on disk.
#!   Such ordering is don on-thy-fly for upto \code{autoindexorder-1} columns that need ordering.
#!   For \code{autoindexorder} o more columns we do the batched ordering upfront with \code{\link{ffindexorder}} and then re-use it in each call to \code{\link{ffindexget}} resp. \code{\link{ffindexset}}.
#! }
#! \value{
#!   Function \code{ffdfindexget} returns a ffdf data.frame with those rows selected by the ff \code{index} vector.
#!   \cr
#!   Function \code{ffdfindexset} returns \code{x} with those rows replaced that had been requested by \code{index} and \code{value}.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{Extract.ff}}, \code{\link{ffindexget}}, \code{\link{ffindexorder}}
#! }
#! \examples{
#! message("ff integer subscripts with ffdf return/assign values")
#! x <- ff(factor(letters))
#! y <- ff(1:26)
#! d <- ffdf(x,y)
#! i <- ff(2:9)
#! di <- d[i,]
#! di
#! d[i,] <- di
#! message("ff integer subscripts: more control with ffindexget/ffindexset")
#! di <- ffdfindexget(d, i, FF_RETURN=di)
#! d <- ffdfindexset(d, i, di)
#! rm(x, y, d, i, di)
#! gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }




ffdfindexget <- function(
  x
, index
, indexorder = NULL
, autoindexorder = 3
, FF_RETURN = NULL
, BATCHSIZE   = NULL
, BATCHBYTES  = getOption("ffmaxbytes")
, VERBOSE     = FALSE
){

  stopifnot(is.ffdf(x))
  stopifnot(is.null(rownames(x)))
  stopifnot(is.ff(index))
  N <- nrow(x)
  n <- length(index)
  stopifnot(vmode(index)=="integer")

  p <- physical(x)
  np <- length(p)
  b <- integer(np)

  v <- character(np)
  dc <- integer(np)

  for (i in 1:np){
    y <- p[[i]]
    if(!is.null(vw(y)))
      stop("subscripting with ff does not support vw")

    v[i] <- vmode(y)
    d <- dim(y)
    if (is.null(d))
      dc[i] <- 0L
    else
      dc[i] <- d[[2]]
  }

  recindbytes <- .rambytes["integer"]
  recvalbytes <- .rambytes[v]
  theobytes.inram <- ifelse(dc, 2*recvalbytes*dc+2*recindbytes, 2*recvalbytes+2*recindbytes)  # verify
  theobytes.chunk <- ifelse(dc, 2*(2*recvalbytes*dc+2*recindbytes), recvalbytes+2*recindbytes)    # verify
  # xx in the case WITH dim AMD chunked AND WITHOUT indexorder we actually use more RAM in unsort.ahi, gc() does not pop-in
  # can we improve on
  # ret <- unsort.ahi(ret, index, ixre)
  # in [.ff_array
  # ?

  if (is.null(BATCHSIZE))
    inram <- n*theobytes.inram<=BATCHBYTES & n>=N
  else
    inram <- n*theobytes.inram<=BATCHBYTES & n>=N & BATCHSIZE>=n

  theobytes <- ifelse(inram, theobytes.inram, theobytes.chunk)

  b <- bbatch(n, as.integer(BATCHBYTES / theobytes))$b

  if (!is.null(BATCHSIZE) && !all(inram))
    b[!inram] <- pmin(BATCHSIZE, b[!inram])

  # if we have many columns we store the chunkwise ordered index (to avoid redundant ordering)
  if (is.null(indexorder)){
    BATCHSIZE <- min(b)
    if ( sum(!inram) >= autoindexorder ){
      indexorder <- ffindexorder(
        index
      , BATCHSIZE
      , VERBOSE     = VERBOSE
      )
      b[!inram] <- BATCHSIZE
    }
  }else{
    stopifnot(is.null(vw(indexorder)))
    stopifnot(vmode(indexorder)=="integer")
    stopifnot(length(indexorder)==n)
    BATCHSIZE <- attr(indexorder, "BATCHSIZE")
    b[!inram] <- BATCHSIZE
  }

  totaltheobytes <- b * theobytes

  if (is.null(FF_RETURN)){
    FF_RETURN <- clone(x, initdata=NULL, nrow=n)
    physical <- physical(FF_RETURN)
  }else{
    physical <- physical(FF_RETURN)
    if(length(physical)!=length(p))
      stop("FF_RETURN not suitable because has different number of physical components")
    for (i in 1:np){
      y <- p[[i]]
      z <- physical[[i]]
      if (!ffsuitable(z, y, FF_ATTR=list(length=n)))
        stop("FF_RETURN not suitable because one physical component is not suitable")
    }
  }

  if (VERBOSE){
    cat(if (is.null(indexorder)) "NOT ", "USING indexorder\n")
    print(data.frame(BATCHES=ceiling(as.double(n)/b), BATCHSIZE=b, InRAM=inram, TheoreticBytes=theobytes, TotalTheoreticBytes=totaltheobytes, DimCols=dc))
  }

  #gc();initram <- memory.size()

  ramindex <- NULL  # initialize since we query it for being is.null()
  for (i in 1:np){
    y <- p[[i]]
    if (dc[i]){  # if object has dimension
      if (inram[i]){
        if(is.null(ramindex))
          ramindex <- read.ff(index, 1L, n)
        physical[[i]][,] <- p[[i]][,,drop=FALSE][ramindex,,drop=FALSE]
        #cat("expected", totaltheobytes[i]/(1024^2), "got"); print(memory.size()-initram)
      }else{
        ramindex <- NULL  # kill RAM reminder
        for (j in chunk(1, n, b[i], method="seq")){
          if (is.null(indexorder)){
            physical[[i]][j,] <- y[index[j],,drop=FALSE]
          }else{
            o <- indexorder[j]+1L
            # the following single-liner fails because of the drop=FALSE symmetry problem
            #physical[[i]][j,,drop=FALSE][o,,drop=FALSE] <- y[index[j][o],,drop=FALSE]
            z <- matrix(do.call(.rammode[vmode(y)], list(1)), nrow=sum(j), ncol=dc[i])
            z[o,] <- y[index[j][o],]
            rm(o)
            physical[[i]][j,] <- z
            rm(z)
          }
          #cat("expected", totaltheobytes[i]/(1024^2), "got");print(memory.size()-initram)
        }
      }


    }else{

      ramindex <- NULL  # kill RAM reminder

      physical[[i]] <- ffindexget(
        y
      , index
      , indexorder = indexorder
      , FF_RETURN = physical[[i]]
      , BATCHSIZE   = b[i]
      , BATCHBYTES  = BATCHBYTES
      , VERBOSE     = FALSE
      )
      #cat("expected", totaltheobytes[i]/(1024^2), "got");print(memory.size()-initram)

    }

  }

  cl <- oldClass(FF_RETURN)
  oldClass(FF_RETURN) <- NULL
  #FF_RETURN$virtual <- virtual
  FF_RETURN$physical <- physical
  oldClass(FF_RETURN) <- cl

  FF_RETURN
}





ffdfindexset <- function(
  x
, index
, value
, indexorder = NULL
, autoindexorder = 3
, BATCHSIZE   = NULL
, BATCHBYTES  = getOption("ffmaxbytes")
, VERBOSE     = FALSE
){

  stopifnot(is.ffdf(x))
  stopifnot(is.ffdf(value))
  stopifnot(is.ff(index))
  stopifnot(vmode(index)=="integer")
  N <- nrow(x)
  n <- length(index)
  stopifnot(nrow(value)==n)

  p <- physical(x)
  physical <- physical(value)
  np <- length(p)
  if(length(physical)!=np)
    stop("value not suitable because has different number of physical components")

  b <- integer(np)
  v <- character(np)
  dc <- integer(np)
  for (i in 1:np){
    y <- p[[i]]
    if(!is.null(vw(y)))
      stop("subscripting with ff does not support vw")
    if (!ffsuitable(physical[[i]], y, FF_ATTR=list(length=n)))
      stop("FF_RETURN not suitable because one physical component is not suitable")

    v[i] <- vmode(y)
    d <- dim(y)
    if (is.null(d))
      dc[i] <- 0L
    else
      dc[i] <- d[[2]]
  }

  recindbytes <- .rambytes["integer"]
  recvalbytes <- .rambytes[v]
  theobytes.inram <- ifelse(dc, 2*recvalbytes*dc+2*recindbytes, 2*recvalbytes+2*recindbytes)  # verify
  theobytes.chunk <- ifelse(dc, 2*(2*recvalbytes*dc+2*recindbytes), recvalbytes+2*recindbytes)    # verify
  # xx in the case WITH dim AMD chunked AND WITHOUT indexorder we actually use more RAM in unsort.ahi, gc() does not pop-in
  # can we improve on
  # ret <- unsort.ahi(ret, index, ixre)
  # in [.ff_array
  # ?

  if (is.null(BATCHSIZE))
    inram <- n*theobytes.inram<=BATCHBYTES & n>=N
  else
    inram <- n*theobytes.inram<=BATCHBYTES & n>=N & BATCHSIZE>=n

  theobytes <- ifelse(inram, theobytes.inram, theobytes.chunk)

  b <- bbatch(n, as.integer(BATCHBYTES / theobytes))$b

  if (!is.null(BATCHSIZE) && !all(inram))
    b[!inram] <- pmin(BATCHSIZE, b[!inram])

  # if we have many columns we store the chunkwise ordered index (to avoid redundant ordering)
  if (is.null(indexorder)){
    BATCHSIZE <- min(b)
    if ( sum(!inram) >= autoindexorder ){
      indexorder <- ffindexorder(
        index
      , BATCHSIZE
      , VERBOSE     = VERBOSE
      )
      b[!inram] <- BATCHSIZE
    }
  }else{
    stopifnot(is.null(vw(indexorder)))
    stopifnot(vmode(indexorder)=="integer")
    stopifnot(length(indexorder)==n)
    BATCHSIZE <- attr(indexorder, "BATCHSIZE")
    b[!inram] <- BATCHSIZE
  }

  totaltheobytes <- b * theobytes

  if (VERBOSE){
    cat(if (is.null(indexorder)) "NOT ", "USING indexorder\n")
    print(data.frame(BATCHES=ceiling(as.double(n)/b), BATCHSIZE=b, InRAM=inram, TheoreticBytes=theobytes, TotalTheoreticBytes=totaltheobytes, DimCols=dc))
  }

  #gc();initram <- memory.size()

  ramindex <- NULL  # initialize since we query it for being is.null()
  for (i in 1:np){
    y <- p[[i]]

    if (dc[i]){  # if object has dimension

      if (inram[i]){
        if(is.null(ramindex))
          ramindex <- read.ff(index, 1L, n)
        p[[i]][,][ramindex,] <- physical[[i]][,,drop=FALSE]
        #cat("expected", totaltheobytes[i]/(1024^2), "got"); print(memory.size()-initram)
      }else{
        ramindex <- NULL  # kill RAM reminder
        for (j in chunk(1, n, b[i], method="seq")){
          if (is.null(indexorder)){
            y[index[j],] <- physical[[i]][j,]
          }else{
            o <- indexorder[j]+1L
            y[index[j][o],] <- physical[[i]][j,,drop=FALSE][o,]
            rm(o)
          }
          #cat("expected", totaltheobytes[i]/(1024^2), "got");print(memory.size()-initram)
        }
      }


    }else{

      ramindex <- NULL  # kill RAM reminder

      # this needs less RAM and is faster compared to a hi method - although we need to read the index data for each column

      ffindexset(
        y
      , index
      , physical[[i]]
      , indexorder = indexorder
      , BATCHSIZE   = b[i]
      , BATCHBYTES  = BATCHBYTES
      , VERBOSE     = FALSE
      )
      #cat("expected", totaltheobytes[i]/(1024^2), "got");print(memory.size()-initram)

    }

  }

  x
}





#! \name{ffdfsort}
#! \alias{dfsort}
#! \alias{dforder}
#! \alias{ramdfsort}
#! \alias{ramdforder}
#! \alias{ffdfsort}
#! \alias{ffdforder}
#! \title{
#! Sorting: convenience wrappers for data.frames
#! }
#! \description{
#! These functions allow convenient sorting and ordering of collections of (ff) vectors organized in (ffdf) data.frames
#! }
#! \usage{
#! dforder(x, ...)
#! dfsort(x, ...)
#! ramdforder(x, ...)
#! ramdfsort(x, ...)
#! ffdforder(x, ...)
#! ffdfsort(x, ...)
#! }
#! \arguments{
#!   \item{x}{
#!   a \code{\link{data.frame}} (for \code{dforder, dfsort, ramorder, ramsort}) or an \code{\link{ffdf}} object (for \code{ffdforder, ffdfsort})
#! }
#!   \item{\dots}{
#!   further arguments passed to \code{\link{sort}}, \code{\link{ramsort}} or \code{\link{ffsort}} (for objects with one column)
#!   or passed to \code{\link{order}}, \code{\link{ramorder}} or \code{\link{fforder}} (for objects with mulitple columns)
#! }
#! }
#! \value{
#!  the order functions return an (ff) vector of integer order positions, the sort functions return a sorted clone of the (ffdf) input data.frame
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#!
#! \seealso{
#!     \code{\link{sort}}, \code{\link{ramsort}} or \code{\link{ffsort}} \cr
#!     \code{\link{order}}, \code{\link{ramorder}} or \code{\link{fforder}}
#! }
#! \examples{
#!    x <- ff(sample(1e5, 1e6, TRUE))
#!    y <- ff(sample(1e5, 1e6, TRUE))
#!    z <- ff(sample(1e5, 1e6, TRUE))
#!    d <- ffdf(x, y, z)
#!    d2 <- ffdfsort(d)
#!    d2
#!    d
#!    d2 <- d[1:2]
#!    i <- ffdforder(d2)
#!    d[i,]
#!    rm(x, y, z, i, d, d2)
#!    gc()
#! }
#! \keyword{univar}
#! \keyword{manip}
#! \keyword{arith}
#! \keyword{ IO }
#! \keyword{ data }




dforder <- function(x, ...){
  k <- ncol(x)
  l <- vector("list", k)
  for (i in 1:k)
    l[[i]] <- x[[i]]
  do.call("order", c(l, ...))
}


dfsort <- function(
  x
, ...
)
{
  if (ncol(x)==1){
    x[[1]] <- sort(x[[1]], ...)
    x
  }else{
    o <- dforder(x
    , ...
    )
    x[o,,drop=FALSE]
  }
}

ramdforder <- function(x, ...){
  K <- ncol(x)
  i <- 1:nrow(x)
  for (k in K:1){
    ramorder(x[[k]], i)
  }
  i
}

ramdfsort <- function(
  x
, ...
)
{
  if (ncol(x)==1){
    x[[1]] <- ramsort(x[[1]], ...)
    x
  }else{
    o <- ramdforder(x, ...)
    x[o,,drop=FALSE]
  }
}



ffdforder <- function(x, ...){
  k <- ncol(x)
  l <- vector("list", k)
  for (i in 1:k)
    l[[i]] <- x[[i]]
  do.call("fforder", c(l, ...))
}

ffdfsort <- function(
  x
, ...
)
{
  if (ncol(x)==1){
    x[[1]] <- ffsort(x[[1]]
    , inplace=FALSE
    , ...
    )
    x
  }else{
    o <- ffdforder(x
    , use.index = FALSE
    , ...
    )
    x[o,,drop=FALSE]
  }
}








#! \name{regtest.fforder}
#! \alias{regtest.fforder}
#! \title{
#! Sorting: regression tests
#! }
#! \description{
#!   Some tests verfying the correctness of the sorting routines
#! }
#! \usage{
#! regtest.fforder(n = 100)
#! }
#! \arguments{
#!   \item{n}{
#!   size of vector to be sorted
#! }
#! }
#! \details{
#!  stops in case of an error
#! }
#! \value{
#!  Invisible()
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{ramsort}}
#! }
#! \examples{
#!   regtest.fforder()
#!
#!  \dontrun{
#!     n <- 5e6
#!     message("performance comparison at n=", n, "")
#!
#!     message("sorting doubles")
#!     x <- y <- as.double(runif(n))
#!
#!     x[] <- y
#!     system.time(sort(x))[3]
#!     x[] <- y
#!     system.time(shellsort(x))[3]
#!     x[] <- y
#!     system.time(shellsort(x, has.na=FALSE))[3]
#!     x[] <- y
#!     system.time(mergesort(x))[3]
#!     x[] <- y
#!     system.time(mergesort(x, has.na=FALSE))[3]
#!
#!     x[] <- y
#!     system.time(sort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE, has.na=FALSE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE, has.na=FALSE))[3]
#!
#!
#!     x <- y <- as.double(sample(c(rep(NA, n/2), runif(n/2))))
#!
#!     x[] <- y
#!     system.time(sort(x))[3]
#!     x[] <- y
#!     system.time(shellsort(x))[3]
#!     x[] <- y
#!     system.time(mergesort(x))[3]
#!
#!     x[] <- y
#!     system.time(sort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE))[3]
#!
#!
#!
#!     x <- y <- sort(as.double(runif(n)))
#!
#!     x[] <- y
#!     system.time(sort(x))  # only here R is faster because R checks for beeing sorted
#!     x[] <- y
#!     system.time(shellsort(x))[3]
#!     x[] <- y
#!     system.time(shellsort(x, has.na=FALSE))[3]
#!     x[] <- y
#!     system.time(mergesort(x))[3]
#!     x[] <- y
#!     system.time(mergesort(x, has.na=FALSE))[3]
#!
#!     x[] <- y
#!     system.time(sort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE, has.na=FALSE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE, has.na=FALSE))[3]
#!
#!     y <- rev(y)
#!     x[] <- y
#!     system.time(sort(x))[3]
#!     x[] <- y
#!     system.time(shellsort(x))[3]
#!     x[] <- y
#!     system.time(shellsort(x, has.na=FALSE))[3]
#!     x[] <- y
#!     system.time(mergesort(x))[3]
#!     x[] <- y
#!     system.time(mergesort(x, has.na=FALSE))[3]
#!
#!     x[] <- y
#!     system.time(sort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(shellsort(x, decreasing=TRUE, has.na=FALSE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE))[3]
#!     x[] <- y
#!     system.time(mergesort(x, decreasing=TRUE, has.na=FALSE))[3]
#!
#!     rm(x,y)
#!
#!
#!     message("ordering doubles")
#!
#!     x <- as.double(runif(n))
#!     system.time(order(x))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i))[3]
#!
#!     x <- as.double(sample(c(rep(NA, n/2), runif(n/2))))
#!     system.time(order(x))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i))[3]
#!
#!     x <- as.double(sort(runif(n)))
#!     system.time(order(x))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i))[3]
#!
#!     x <- rev(x)
#!     system.time(order(x))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i))[3]
#!
#!
#!     x <- as.double(runif(n))
#!     system.time(order(x, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i, decreasing=TRUE))[3]
#!
#!     x <- as.double(sample(c(rep(NA, n/2), runif(n/2))))
#!     system.time(order(x, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i, decreasing=TRUE))[3]
#!
#!     x <- as.double(sort(runif(n)))
#!     system.time(order(x, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i, decreasing=TRUE))[3]
#!
#!     x <- rev(x)
#!     system.time(order(x, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE))[3]
#!     i <- 1:n
#!     system.time(shellorder(x, i, decreasing=TRUE, stabilize=TRUE))[3]
#!     i <- 1:n
#!     system.time(mergeorder(x, i, decreasing=TRUE))[3]
#!
#!
#!     keys <- c("short","ushort")
#!     for (v in c("integer", keys)){
#!
#!       if (v \%in\% keys){
#!         k <- .vmax[v]-.vmin[v]+1L
#!         if (is.na(.vNA[v])){
#!           y <- sample(c(rep(NA, k), .vmin[v]:.vmax[v]), n, TRUE)
#!         }else{
#!           y <- sample(.vmin[v]:.vmax[v], n, TRUE)
#!         }
#!       }else{
#!         k <- .Machine$integer.max
#!         y <- sample(k, n, TRUE)
#!       }
#!
#!       message("sorting ",v)
#!
#!       x <- y
#!       message("sort(x) ", system.time(sort(x))[3])
#!       x <- y
#!       message("shellsort(x) ", system.time(shellsort(x))[3])
#!       x <- y
#!       message("mergesort(x) ", system.time(mergesort(x))[3])
#!       x <- y
#!       message("radixsort(x) ", system.time(radixsort(x))[3])
#!       if (v \%in\% keys){
#!         x <- y
#!         message("keysort(x) ", system.time(keysort(x))[3])
#!         x <- y
#!         message("keysort(x, keyrange=c(.vmin[v],.vmax[v])) "
#! , system.time(keysort(x, keyrange=c(.vmin[v],.vmax[v])))[3])
#!       }
#!
#!       if (!is.na(.vNA[v])){
#!         x <- y
#!         message("shellsort(x, has.na=FALSE) ", system.time(shellsort(x, has.na=FALSE))[3])
#!         x <- y
#!         message("mergesort(x, has.na=FALSE) ", system.time(mergesort(x, has.na=FALSE))[3])
#!         x <- y
#!         message("radixsort(x, has.na=FALSE) ", system.time(radixsort(x, has.na=FALSE))[3])
#!         if (v \%in\% keys){
#!           x <- y
#!           message("keysort(x, has.na=FALSE) ", system.time(keysort(x, has.na=FALSE))[3])
#!           x <- y
#!           message("keysort(x, has.na=FALSE, keyrange=c(.vmin[v],.vmax[v])) "
#! , system.time(keysort(x, has.na=FALSE, keyrange=c(.vmin[v],.vmax[v])))[3])
#!         }
#!       }
#!
#!
#!       message("ordering",v)
#!
#!       x[] <- y
#!       i <- 1:n
#!       message("order(x) ", system.time(order(x))[3])
#!       x[] <- y
#!       i <- 1:n
#!       message("shellorder(x, i) ", system.time(shellorder(x, i))[3])
#!       x[] <- y
#!       i <- 1:n
#!       message("mergeorder(x, i) ", system.time(mergeorder(x, i))[3])
#!       x[] <- y
#!       i <- 1:n
#!       message("radixorder(x, i) ", system.time(radixorder(x, i))[3])
#!       if (v \%in\% keys){
#!         x[] <- y
#!         i <- 1:n
#!         message("keyorder(x, i) ", system.time(keyorder(x, i))[3])
#!         x[] <- y
#!         i <- 1:n
#!         message("keyorder(x, i, keyrange=c(.vmin[v],.vmax[v])) "
#! , system.time(keyorder(x, i, keyrange=c(.vmin[v],.vmax[v])))[3])
#!       }
#!
#!       if (!is.na(.vNA[v])){
#!         x[] <- y
#!         i <- 1:n
#!         message("shellorder(x, i, has.na=FALSE) ", system.time(shellorder(x, i, has.na=FALSE))[3])
#!         x[] <- y
#!         i <- 1:n
#!         message("mergeorder(x, i, has.na=FALSE) ", system.time(mergeorder(x, i, has.na=FALSE))[3])
#!         x[] <- y
#!         i <- 1:n
#!         message("radixorder(x, i, has.na=FALSE) ", system.time(radixorder(x, i, has.na=FALSE))[3])
#!         if (v \%in\% keys){
#!           x[] <- y
#!           i <- 1:n
#!           message("keyorder(x, i, has.na=FALSE) ", system.time(keyorder(x, i, has.na=FALSE))[3])
#!           x[] <- y
#!           i <- 1:n
#!           message("keyorder(x, i, has.na=FALSE, keyrange=c(.vmin[v],.vmax[v])) "
#! , system.time(keyorder(x, i, has.na=FALSE, keyrange=c(.vmin[v],.vmax[v])))[3])
#!         }
#!       }
#!
#!     }
#!   }
#! }
#! \keyword{univar}
#! \keyword{manip}
#! \keyword{arith}


regtest.fforder <- function(n=100){

	# insituorder integer

	#x <- as.integer(runif(n, 1, n))
	#o <- order(x)
	#insituorder(x, o)
	#stopifnot(identical(x, sort(x)))

	# applyorder double

	#x <- runif(n, 1, n)
	#o <- order(x)
	#insituorder(x, o)
	#stopifnot(identical(x, sort(x)))

	# shellsort integer

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x, has.na=FALSE)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x, decreasing=TRUE, has.na=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x)))
	stopifnot(identical(x, sort(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x, na.last=FALSE)))
	stopifnot(identical(x, sort(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=FALSE)))


	# shellsort double

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x)))
	stopifnot(identical(x, sort(x)))

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x, has.na=FALSE)))
	stopifnot(identical(x, sort(x)))

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, shellsort(x, decreasing=TRUE, has.na=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x)))
	stopifnot(identical(x, sort(x, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x, na.last=FALSE)))
	stopifnot(identical(x, sort(x, na.last=FALSE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), shellsort(x, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=FALSE)))





	# mergesort integer

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x, has.na=FALSE)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x, decreasing=TRUE, has.na=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x)))
	stopifnot(identical(x, sort(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x, na.last=FALSE)))
	stopifnot(identical(x, sort(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=FALSE)))


	# mergesort double

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x)))
	stopifnot(identical(x, sort(x)))

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x, has.na=FALSE)))
	stopifnot(identical(x, sort(x)))

	x <- as.double(runif(n, 1, n))
	stopifnot(identical(0L, mergesort(x, decreasing=TRUE, has.na=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x)))
	stopifnot(identical(x, sort(x, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x, na.last=FALSE)))
	stopifnot(identical(x, sort(x, na.last=FALSE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), mergesort(x, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=FALSE)))







	# shellorder integer without stabilize

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i)))
	stopifnot(identical(x[i], sort(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE)))
	stopifnot(identical(x[i], sort(x,decreasing=TRUE)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,has.na=FALSE)))
	stopifnot(identical(x[i], sort(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE,has.na=FALSE)))
	stopifnot(identical(x[i], sort(x,decreasing=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i)))
	stopifnot(identical(x[i], sort(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE)))
	stopifnot(identical(x[i], sort(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, na.last=FALSE)))
	stopifnot(identical(x[i], sort(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x[i], sort(x, decreasing=TRUE, na.last=FALSE)))


	# shellorder integer with stabilize

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,stabilize=TRUE)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE,stabilize=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,has.na=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE,has.na=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i,stabilize=TRUE)))
	stopifnot(identical(i, order(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE,stabilize=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, na.last=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE, na.last=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=FALSE)))




	# shellorder double without stabilize

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i)))
	stopifnot(identical(x[i], sort(x)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE)))
	stopifnot(identical(x[i], sort(x,decreasing=TRUE)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,has.na=FALSE)))
	stopifnot(identical(x[i], sort(x)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE,has.na=FALSE)))
	stopifnot(identical(x[i], sort(x,decreasing=TRUE)))


	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i)))
	stopifnot(identical(x[i], sort(x, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE)))
	stopifnot(identical(x[i], sort(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, na.last=FALSE)))
	stopifnot(identical(x[i], sort(x, na.last=FALSE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x[i], sort(x, decreasing=TRUE, na.last=FALSE)))


	# shellorder double with stabilize

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,stabilize=TRUE)))
	stopifnot(identical(i, order(x)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE,stabilize=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,has.na=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, shellorder(x,i,decreasing=TRUE,has.na=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))


	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i,stabilize=TRUE)))
	stopifnot(identical(i, order(x, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE,stabilize=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, na.last=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x, na.last=FALSE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), shellorder(x,i, decreasing=TRUE, na.last=FALSE,stabilize=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=FALSE)))



	# mergeorder integer

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i,decreasing=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i,has.na=FALSE)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i,decreasing=TRUE,has.na=FALSE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), nNA <- mergeorder(x,i)))
	stopifnot(identical(i, order(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), mergeorder(x,i, decreasing=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), mergeorder(x,i, na.last=FALSE)))
	stopifnot(identical(i, order(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), mergeorder(x,i, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=FALSE)))


	# mergeorder double

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i)))
	stopifnot(identical(i, order(x)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i,decreasing=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i,has.na=FALSE)))
	stopifnot(identical(i, order(x)))

	x <- as.double(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, mergeorder(x,i,decreasing=TRUE,has.na=FALSE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))


	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), nNA <- mergeorder(x,i)))
	stopifnot(identical(i, order(x, na.last=TRUE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), mergeorder(x,i, decreasing=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), mergeorder(x,i, na.last=FALSE)))
	stopifnot(identical(i, order(x, na.last=FALSE)))

	x <- as.double(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), mergeorder(x,i, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=FALSE)))


	# keysort integer

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	stopifnot(identical(0L, keysort(x)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	stopifnot(identical(0L, keysort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	stopifnot(identical(0L, keysort(x, has.na=FALSE)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	stopifnot(identical(0L, keysort(x, decreasing=TRUE, has.na=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), keysort(x)))
	stopifnot(identical(x, sort(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), keysort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), keysort(x, na.last=FALSE)))
	stopifnot(identical(x, sort(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), keysort(x, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=FALSE)))


	# keyorder integer

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	i <- 1:n
	stopifnot(identical(0L, keyorder(x,i)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	i <- 1:n
	stopifnot(identical(0L, keyorder(x,i,decreasing=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	i <- 1:n
	stopifnot(identical(0L, keyorder(x,i,has.na=FALSE)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, .vmin["short"], .vmax["short"]))
	i <- 1:n
	stopifnot(identical(0L, keyorder(x,i,decreasing=TRUE,has.na=FALSE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), nNA <- keyorder(x,i)))
	stopifnot(identical(i, order(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), keyorder(x,i, decreasing=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), keyorder(x,i, na.last=FALSE)))
	stopifnot(identical(i, order(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["short"], .vmax["short"])), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), keyorder(x,i, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=FALSE)))




	# radixsort integer

	x <- as.integer(runif(n, .vmin["integer"], .vmax["integer"]))
	stopifnot(identical(0L, radixsort(x)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, .vmin["integer"], .vmax["integer"]))
	stopifnot(identical(0L, radixsort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(runif(n, .vmin["integer"], .vmax["integer"]))
	stopifnot(identical(0L, radixsort(x, has.na=FALSE)))
	stopifnot(identical(x, sort(x)))

	x <- as.integer(runif(n, .vmin["integer"], .vmax["integer"]))
	stopifnot(identical(0L, radixsort(x, decreasing=TRUE, has.na=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["integer"], .vmax["integer"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), radixsort(x)))
	stopifnot(identical(x, sort(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["integer"], .vmax["integer"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), radixsort(x, decreasing=TRUE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["integer"], .vmax["integer"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), radixsort(x, na.last=FALSE)))
	stopifnot(identical(x, sort(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, .vmin["integer"], .vmax["integer"])), 2*n, TRUE))
	stopifnot(identical(sum(is.na(x)), radixsort(x, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(x, sort(x, decreasing=TRUE, na.last=FALSE)))



	# radixorder integer

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, radixorder(x,i)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, radixorder(x,i,decreasing=TRUE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, radixorder(x,i,has.na=FALSE)))
	stopifnot(identical(i, order(x)))

	x <- as.integer(runif(n, 1, n))
	i <- 1:n
	stopifnot(identical(0L, radixorder(x,i,decreasing=TRUE,has.na=FALSE)))
	stopifnot(identical(i, order(x,decreasing=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), nNA <- radixorder(x,i)))
	stopifnot(identical(i, order(x, na.last=TRUE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), radixorder(x,i, decreasing=TRUE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=TRUE)))


	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), radixorder(x,i, na.last=FALSE)))
	stopifnot(identical(i, order(x, na.last=FALSE)))

	x <- as.integer(sample(c(rep(NA, n), runif(n, 1, n)), 2*n, TRUE))
	i <- 1:(2*n)
	stopifnot(identical(sum(is.na(x)), radixorder(x,i, decreasing=TRUE, na.last=FALSE)))
	stopifnot(identical(i, order(x, decreasing=TRUE, na.last=FALSE)))



	n <- 1e2

	x <- sample(c(Inf, Inf, 5,5, 0,0, NA,NA, NaN, NaN, -5,-5, -Inf,-Inf), n, TRUE)
	i <- 1:length(x)
	mergeorder(x, i)
	stopifnot(identical(order(x), i))

	x <- sample(c(Inf, Inf, 5,5, 0,0, NA,NA, NaN, NaN, -5,-5, -Inf,-Inf), n, TRUE)
	i <- 1:length(x)
	shellorder(x, i, stabilize=TRUE)
	stopifnot(identical(order(x), i))

	x <- sample(c(Inf, Inf, 5,5, 0,0, NA,NA, -5,-5, -Inf,-Inf), n, TRUE)
	i <- 1:length(x)
	shellorder(x, i)
	# since shellorder is not stable, the result is NOT identical with order()
	stopifnot(identical(sort(x, na.last=TRUE), x[i]))



	x <- as.integer(sample(c(5,5, 0,0, NA,NA, -5,-5), n, TRUE))
	i <- 1:length(x)
	shellorder(x, i)
	# since shellorder is not stable, the result is NOT identical with order()
	stopifnot(identical(sort(x, na.last=TRUE), x[i]))

	x <- as.integer(sample(c(5,5, 0,0, NA,NA, -5,-5), n, TRUE))
	i <- 1:length(x)
	shellorder(x, i, stabilize=TRUE)
	stopifnot(identical(order(x), i))

	x <- as.integer(sample(c(5,5, 0,0, NA,NA, -5,-5), n, TRUE))
	i <- 1:length(x)
	mergeorder(x, i)
	stopifnot(identical(order(x), i))

	x <- as.integer(sample(c(5,5, 0,0, NA,NA, -5,-5), n, TRUE))
	i <- 1:length(x)
	radixorder(x, i)
	stopifnot(identical(order(x), i))

	x <- as.integer(sample(c(5,5, 0,0, NA,NA, -5,-5), n, TRUE))
	i <- 1:length(x)
	keyorder(x, i)
stopifnot(identical(order(x), i))


}





if (FALSE){



	# test ffsort

	library(ff)
	n <- 26e4L
	ordersize <- 26e4L
	has.na <- TRUE
	na.last <- TRUE
	decreasing <- FALSE

	x <- ff(factor(rev(letters)), vmode="short", length=n)
	system.time(ffsort(x, verbose=TRUE))

	x <- ff(26:1, vmode="short", length=n)
	system.time(ffsort(x, verbose=TRUE))

	x <- ff(26:1, vmode="integer", length=n)
	system.time(ffsort(x, verbose=TRUE))









	library(ff)
	n <- 52e6L
	has.na <- TRUE
	na.last <- TRUE
	decreasing <- FALSE

	i <- ff(vmode="integer", length=n)

	x <- ff(26:1, vmode="short", length=n)
	y <- clone(x)
	for (j in chunk(i))i[j] <- j[[1]]:j[[2]]
	system.time(fforder(x, i, verbose=TRUE))
	j <- 1:n; system.time(ramorder(y[], j))
	identical(i[], j)

	x <- ff(26:1, vmode="integer", length=n)
	y <- clone(x)
	for (j in chunk(i))i[j] <- j[[1]]:j[[2]]
	system.time(fforder(x, i, verbose=TRUE))
	j <- 1:n; system.time(ramorder(y[], j))
	identical(i[], j)
	system.time(j <- order(y[]))
	identical(i[], j)

	x <- ff(26:1, vmode="double", length=n)
	y <- clone(x)
	for (j in chunk(i))i[j] <- j[[1]]:j[[2]]
	system.time(fforder(x, i, verbose=TRUE))
	j <- 1:n; system.time(ramorder(y[], j))
	identical(i[], j)
	system.time(j <- order(y[]))
	identical(i[], j)

	x <- ff(26:1, vmode="single", length=n)
	y <- clone(x)
	for (j in chunk(i))i[j] <- j[[1]]:j[[2]]
	system.time(fforder(x, i, verbose=TRUE))
	j <- 1:n; system.time(ramorder(y[], j))
	identical(i[], j)

	x <- ff(factor(rev(letters)), length=n)
	y <- clone(x)
	for (j in chunk(i))i[j] <- j[[1]]:j[[2]]
	system.time(fforder(x, i, verbose=TRUE))
	j <- 1:n; system.time(ramorder(y[], j))
	identical(i[], j)
	system.time(j <- order(y[]))
	identical(i[], j)


	library(ff)
	x1 <- sample(5, 25, TRUE)
	x2 <- sample(5, 25, TRUE)
	f1 <- ff(x1)
	f2 <- ff(x2)
	fo <- ff(1:25)

	o <- 1:25
	o1 <- f1[]
	o2 <- f2[]
	ramorder(o2, o)
	ramorder(o1, o)
	cbind(x1, x2)[o,]

	fforder(f2, fo)
	fforder(f1, fo, use.index=TRUE)
	cbind(x1,x2)[fo[],]










	nNA <- .Call("ffkeysort"
				, ffmode_       = .ffmode[vmode(x)]
				, ff_           = attr(x, "physical")
				, left_         = 1L
				, right_        = length(x)
				, keyrange_     = as.integer(keyrange)
				, ordersize_    = as.integer(ordersize)
				, has_na_       = as.logical(has.na)
				, na_last_      = as.logical(na.last)
				, decreasing_   = as.logical(decreasing)
				, PACKAGE="ff"
				)

	# stürzt nicht ab
	x <- ff(sample(as.factor(letters)), length=100000)
	ffsort(x)

	# stürzt ab
	x <- ff(sample(as.factor(letters)), length=100000, vmode="byte")
	ffsort(x)

	x <- ff(sample(as.factor(letters)), length=100000, vmode="ubyte")
	ffsort(x)







	library(ff)

	n <- 1e6L
	o <- 1e6L
	m <- as.integer(32768/2) #1e6
	v <- "double"
	keyrange <- c(1L, 1000L)

	ffx <- ff(vmode=v, length=n, caching="mmeachflush")
	left <- 1L
	right <- n

	s <- 1L
	#for (s in 1:1000){
		cat("seed",s,"\n")
		#sink("d:/tmp/t.txt")
		set.seed(s)

		#for (i in chunk(ffx))ffx[i] <- sample(c(Inf, 1, 0, NA, NaN, -1, -Inf, runif(7)), sum(i), TRUE)
		#for (i in chunk(ffx))ffx[i] <- sample(c(Inf, 1, 0, -1, -Inf, runif(7)), sum(i), TRUE)
		#for (i in chunk(ffx))ffx[i] <- sample(c(keyrange[1]:keyrange[2], rep(NA, diff(keyrange))), sum(i), TRUE)
		#for (i in chunk(ffx))ffx[i] <- sample(keyrange[1]:keyrange[2], sum(i), TRUE)
		#system.time(for (i in chunk(ffx))ffx[i] <- as.integer(runif(sum(i), max=.Machine$integer.max)))
		system.time(for (i in chunk(ffx))ffx[i] <- runif(sum(i), max=.Machine$integer.max))

		system.time(ffx2 <- clone(ffx))
		system.time(
		nNA <- .Call("ffkeysort"
		, ff_           = attr(ffx2, "physical")
		, left_         = left
		, right_        = right
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, has_na_       = TRUE
		, na_last_      = TRUE
		, decreasing_   = TRUE
		, PACKAGE="ff"
		)
		)






		x <- ffx[]
		sum(is.na(x))
		nNA

		system.time(i <- sort(x, decreasing=TRUE, na.last=TRUE))
		identical(ffx2[], i)



		ffy <- ff(vmode=vmode(ffx), length=length(ffx), caching="mmeachflush")

		method <- 0L
		system.time(ffx2 <- clone(ffx))

		system.time(
		.Call("ffsortmerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, auxff_        = attr(ffy, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)

		x <- ffx[]
		sum(is.na(x))
		nNA

		system.time(i <- sort(x, decreasing=TRUE, na.last=FALSE))
		identical(ffx2[], i)


		method <- 1L
		system.time(ffx2 <- clone(ffx))

		system.time(
		.Call("ffsortmerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, auxff_        = attr(ffy, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)

		x <- ffx[]
		sum(is.na(x))
		nNA

		system.time(i <- sort(x, decreasing=TRUE, na.last=FALSE))
		identical(ffx2[], i)


		method <- 2L
		system.time(ffx2 <- clone(ffx))

		system.time(
		.Call("ffsortmerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, auxff_        = attr(ffy, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)

		x <- ffx[]
		sum(is.na(x))
		nNA

		system.time(i <- sort(x, decreasing=TRUE, na.last=FALSE))
		identical(ffx2[], i)


		method <- 3L
		system.time(ffx2 <- clone(ffx))

		system.time(
		.Call("ffsortmerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, auxff_        = attr(ffy, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)

		x <- ffx[]
		sum(is.na(x))
		nNA

		system.time(i <- sort(x, decreasing=TRUE, na.last=FALSE))
		identical(ffx2[], i)








		ffo <- ff(vmode="integer", length=n, caching="mmeachflush")
		ffp <- ff(vmode=vmode(ffo), length=length(ffo), caching="mmeachflush")


		method <- 0L
		for (i in chunk(ffo))ffo[i] <- (i[[1]]):(i[[2]])
		system.time(ffx2 <- clone(ffx))

		system.time(
		nNA <- .Call("ffordermerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, index_        = attr(ffo, "physical")
		, auxff_        = attr(ffy, "physical")
		, auxindex_     = attr(ffp, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_     = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, orderindex_   = TRUE
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)

		#ffx[ffo[]]
		x <- ffx[]
		sum(is.na(x))
		nNA

		#system.time(i <- order(x, decreasing=TRUE, na.last=FALSE))
		system.time({i <- 1:n;mergeorder(x, i, decreasing=TRUE, na.last=FALSE)})
		identical(ffo[], i)



		method <- 1L
		system.time(for (i in chunk(ffo))ffo[i] <- (i[[1]]):(i[[2]]))
		system.time(ffx2 <- clone(ffx))

		system.time(
		nNA <- .Call("ffordermerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, index_        = attr(ffo, "physical")
		, auxff_        = attr(ffy, "physical")
		, auxindex_     = attr(ffp, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, orderindex_   = TRUE
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)

		#ffx[ffo[]]
		x <- ffx[]
		sum(is.na(x))
		nNA

		#system.time(i <- order(x, decreasing=TRUE, na.last=FALSE))
		system.time({i <- 1:n;mergeorder(x, i, decreasing=TRUE, na.last=FALSE)})
		identical(ffo[], i)



		method <- 2L
		for (i in chunk(ffo))ffo[i] <- (i[[1]]):(i[[2]])
		system.time(ffx2 <- clone(ffx))

		system.time(
		nNA <- .Call("ffordermerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, index_        = attr(ffo, "physical")
		, auxff_        = attr(ffy, "physical")
		, auxindex_     = attr(ffp, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, orderindex_   = TRUE
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)
		#sink()

		#ffx[ffo[]]
		x <- ffx[]
		sum(is.na(x))
		nNA

		#system.time(i <- order(x, decreasing=TRUE, na.last=FALSE))
		system.time({i <- 1:n;mergeorder(x, i, decreasing=TRUE, na.last=FALSE)})
		identical(ffo[], i)



		method <- 3L
		for (i in chunk(ffo))ffo[i] <- (i[[1]]):(i[[2]])
		system.time(ffx2 <- clone(ffx))

		system.time(
		nNA <- .Call("ffordermerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, index_        = attr(ffo, "physical")
		, auxff_        = attr(ffy, "physical")
		, auxindex_     = attr(ffp, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = as.integer(keyrange)
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, orderindex_   = TRUE
		, has_na_       = TRUE
		, na_last_      = FALSE #TRUE
		, decreasing_   = TRUE #FALSE
		, PACKAGE="ff"
		)
		)
		#sink()

		#ffx[ffo[]]
		x <- ffx[]
		sum(is.na(x))
		nNA

		#system.time(i <- order(x, decreasing=TRUE, na.last=FALSE))
		system.time({i <- 1:n;mergeorder(x, i, decreasing=TRUE, na.last=FALSE)})
		identical(ffo[], i)







		library(ff)
		n <- 1e8L
		o <- 5e7L
		m <- as.integer(32768/2) #1e6
		v <- "double"
		caching <- "mmeachflush"

		ffx <- ff(vmode=v, length=n, caching=caching)
		left <- 1L
		right <- n

		s <- 1L
		cat("seed",s,"\n")
		set.seed(s)

		#for (i in chunk(ffx))ffx[i] <- sample(c(Inf, 1, 0, NA, NaN, -1, -Inf, runif(7)), sum(i), TRUE)
		#system.time(for (i in chunk(ffx))ffx[i] <- as.integer(runif(sum(i), max=n)))
		system.time(for (i in chunk(ffx))ffx[i] <- runif(sum(i), max=n))

		ffy <- ff(vmode=vmode(ffx), length=length(ffx), caching=caching)
		ffo <- ff(vmode="integer", length=n, caching=caching)
		ffp <- clone(ffo)

		method <- 0L
		system.time(ffx2 <- clone(ffx))
		for (i in chunk(ffo))ffo[i] <- (i[[1]]):(i[[2]])

		system.time(
		nNA <- .Call("ffordermerge"
		, ffmode_       = .ffmode[vmode(ffx2)]
		, ff_           = attr(ffx2, "physical")
		, index_        = attr(ffo, "physical")
		, auxff_        = attr(ffy, "physical")
		, auxindex_     = attr(ffp, "physical")
		, left_         = left
		, right_        = right
		, method_       = as.integer(method)
		, keyrange_       = c(1L, .Machine$integer.max) # NULL
		, ordersize_    = as.integer(o)
		, mergesize_    = as.integer(m)
		, orderindex_   = TRUE
		, has_na_       = TRUE
		, na_last_      = TRUE
		, decreasing_   = FALSE
		, PACKAGE="ff"
		)
		)

		#x <- ffx[]
		#sum(is.na(x))
		#nNA
		#system.time(i <- sort(x, decreasing=FALSE, na.last=TRUE))
		#identical(ffx2[], i)

		#m <- matrix(ffo[], nrow=o)
		#for (i in 1:ncol(m))
		#m[,i] <- sort(m[,i])

		#io <- (0:(n %/% o))*o
		#io <- rep(io, rep(o, length(io)))[1:n]
		#io <- io + 1
		#x <- ffo[]
		#system.time(
		#sort.int(x, index.return=TRUE, method="quick")
		#)
		#x <- ffo[]
		#i <- 1:n
		#system.time({
		#radixorder(x, i, has.na=FALSE)
		#})

		method2 <- 4L

		system.time(
		.Call("ffchunkorder"
		, index_        = attr(ffo, "physical")
		, outindex_        = attr(ffp, "physical")
		, indexsize_  = as.integer(n)
		, method_       = method2  #/* 0=mergeorder 1=shellorder 2=radixorder 4=quickorder */
		, ordersize_    = as.integer(o) #/* int no of elements to be ordered in RAM (must be same as in r_ff_index_chunkorder) */
		)
		)
		#identical(matrix(ffo[ffp[]+io], nrow=o), m)


		#ffy[] <- 0
		system.time(
		.Call("ffindexget"
		, ffmode_       = .ffmode[vmode(ffx)]
		, ffin_         = attr(ffx, "physical")
		, ffout_        = attr(ffy, "physical")
		, index_        = attr(ffo, "physical")
		, auxindex_     = NULL
		, offset_       = 1  # 1 for R2C
		, left_         = 1
		, right_        = n
		, method_       = method2  #/* 0=mergeorder 1=shellorder */
		, ordersize_    = as.integer(o) #/* int no of elements to be ordered in RAM (must be same as in r_ff_index_chunkorder) */
		)
		)
		#identical(ffx2[],ffy[])


		#sink("d:/tmp/t.txt")

		#ffy[] <- 0
		system.time(
		.Call("ffindexget"
		, ffmode_       = .ffmode[vmode(ffx)]
		, ffin_         = attr(ffx, "physical")
		, ffout_        = attr(ffy, "physical")
		, index_        = attr(ffo, "physical")
		, auxindex_     = attr(ffp, "physical")
		, offset_       = 1  # 1 for R2C
		, left_         = 1
		, right_        = n
		, method_       = method2  #/* 0=mergeorder 1=shellorder */
		, ordersize_    = as.integer(o) #/* int no of elements to be ordered in RAM (must be same as in r_ff_index_chunkorder) */
		)
		)
		#sink()
		#identical(ffx2[],ffy[])


		#ffy[] <- 0
		system.time(
		.Call("ffindexget"
		, ffmode_       = .ffmode[vmode(ffx)]
		, ffin_         = attr(ffx, "physical")
		, ffout_        = attr(ffy, "physical")
		, index_        = attr(ffo, "physical")
		, auxindex_     = FALSE
		, offset_       = 1  # 1 for R2C
		, left_         = 1
		, right_        = n
		, method_       = method2  #/* 0=mergeorder 1=shellorder */
		, ordersize_    = as.integer(o) #/* int no of elements to be ordered in RAM (must be same as in r_ff_index_chunkorder) */
		)
		)
		#identical(ffx2[],ffy[])




		#data test(keep=x);
		#   seed = 1;
		#   do i = 1 to 1e8;
		#      call ranuni(seed, x);
		#      output;
		#   end;
		#run;

		#proc sort data=test;
		#by x;
		#run;



		#data test(keep=x y);
		#   seed = 1;
		#   do i = 1 to 1e8;
		#      call ranuni(seed, x);
		#      call ranuni(seed, y);
		#      output;
		#   end;
		#run;

		#proc sort data=test;
		#by x;
		#run;



		#data test(keep=x y z);
		#   seed = 1;
		#   do i = 1 to 1e8;
		#      call ranuni(seed, x);
		#      call ranuni(seed, y);
		#      call ranuni(seed, z);
		#      output;
		#   end;
		#run;

		#proc sort data=test;
		#by x;
		#run;




	#}


		dforder <- function(x, ...){
			k <- ncol(x)
			l <- vector("list", k)
			for (i in 1:k)
				l[[i]] <- x[[i]]
			do.call("order", c(l, ...))
		}


		library(ff)

		N <- c(3e6L, 9e6L, 27e6, 81e6)
		M <- c(2L, 4L, 6L)
		K <- c(1L, 2L)
		tim <- array(NA, dim=c(2, 4, length(N), length(M), length(K)), dimnames=list(c("R","ff")
		, c("o <- order(x)", "xo <- x[o]", "x2 <- clone(xo)", "x2[o] <- xo"), N, M, K))

		ki <- 1L
		mi <- 1L
		ni <- 1L

		for (ki in 1:length(K)){
			for (mi in 1:length(M)){
				for (ni in 1:length(N)){
					k <- K[ki]
					m <- M[mi]
					n <- N[ni]

					vmode <- rep("integer", m)
					l <- list()
					for (j in 1:m){
						l[[j]] <- ff(vmode=vmode[j], length=n)
					}
					for (i in chunk(l[[2]])){
						h <- as.hi(i)
						for (j in 1:m){
							if (j==1)
								x <- sample(1000, sum(i), TRUE)
							else if (j==2)
								x <- runif(sum(i), 1, n)
							l[[j]][h] <- x + j
						}
						rm(x)
					}
					names(l) <- paste("c", 1:m, sep="")
					fd <- do.call("ffdf", l)

					if (TRUE){
						gc()
						tim["ff", 1,ni,mi,ki] <- system.time( fo <- ffdforder(fd[1:k])          )[3]
						tim["ff", 2,ni,mi,ki] <- system.time( fdo <- fd[fo,]                    )[3]
						tim["ff", 3,ni,mi,ki] <- system.time( fd2 <- clone(fdo)                 )[3]
						tim["ff", 4,ni,mi,ki] <- system.time( fd2[fo,] <- fdo                   )[3]
					}
					if (n<=9e6){  # limit to cases where we do not exhaust RAM
						d <- fd[,]
						gc()
						tim["R", 1,ni,mi,ki] <- system.time( o <- dforder(d[1:k]) )[3]
						tim["R", 2,ni,mi,ki] <- system.time( do <- d[o,]          )[3]
						rm(d)
						gc()
						tim["R", 3,ni,mi,ki] <- system.time( d2 <- do             )[3]
						tim["R", 4,ni,mi,ki] <- system.time( d2[o,] <- do         )[3]
						rm(o, do, d2)
					}

					print(tim)
					gc()
				}
			}
		}


}

