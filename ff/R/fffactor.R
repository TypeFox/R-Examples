# some factor support for ff and ffdf objects
# (c) 2009 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2009-09-09
# Last changed: 2009-09-09

# source("d:/mwp/eanalysis/ff/R/fffactor.R")


#! \name{sortLevels}
#! \Rdversion{1.1}
#! \alias{appendLevels}
#! \alias{recodeLevels}
#! \alias{recodeLevels.factor}
#! \alias{recodeLevels.ff}
#! \alias{sortLevels}
#! \alias{sortLevels.factor}
#! \alias{sortLevels.ff}
#! \alias{sortLevels.ffdf}
#! \title{
#!   Factor level manipulation
#! }
#! \description{
#!   \code{appendLevels} combines \code{\link[=levels.ff]{levels}} without sorting such that levels of the first argument will not require re-coding.
#!   \cr
#!   \code{recodeLevels} is a generic for recoding a factor to a desired set of levels - also has a method for large \code{\link{ff}} objects
#!   \cr
#!   \code{sortLevels} is a generic for level sorting and recoding of single factors or of all factors of a \code{\link{ffdf}} dataframe.
#! }
#! \usage{
#! appendLevels(\dots)
#! recodeLevels(x, lev)
#! \method{recodeLevels}{factor}(x, lev)
#! \method{recodeLevels}{ff}(x, lev)
#! sortLevels(x)
#! \method{sortLevels}{factor}(x)
#! \method{sortLevels}{ff}(x)
#! \method{sortLevels}{ffdf}(x)
#! }
#! \arguments{
#!   \item{\dots}{character vector of levels or \code{\link[=is.factor.ff]{is.factor}} objects from which the level attribute is taken}
#!   \item{x}{a \code{\link{factor}} or \code{\link{ff}} factor or a \code{\link{ffdf}} dataframe (\code{sortLevels} only) }
#!   \item{lev}{a character vector of levels}
#! }
#! \details{
#!   When reading a long file with categorical columns the final set of factor levels is only known once the complete file has been read.
#!   When a file is so large that we read it in chunks, the new levels need to be added incrementally.
#!   \code{\link[base]{rbind.data.frame}} sorts combined levels, which requires recoding. For \code{\link{ff}} factors this would require recoding of all previous chunks at the next chunk - potentially on disk, which is too expensive.
#!   Therefore \code{\link{read.table.ffdf}} will simply \code{appendLevels} without sorting, and the \code{recodeLevels} and \code{sortLevels} generics provide a convenient means for sorting and recoding levels after all chunks have been read.
#! }
#! \value{
#!   \code{appendLevels} returns a vector of combined levels, \code{recodeLevels} and \code{sortLevels} return the input object with changed levels. Do read the note!
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \note{
#!   You need to re-assign the return value not only for ram- but also for ff-objects. Remember ff's hybrid copying semantics: \code{\link{LimWarn}}.
#!   If you forget to re-assign the returned object, you will end up with ff objects that have their integer codes re-coded to the new levels but still carry the old levels as a \code{\link[=virtual.ff]{virtual}} attribute.
#! }
#! \seealso{
#!   \code{\link{read.table.ffdf}}, \code{\link{levels.ff}}
#! }
#! \examples{
#!   message("Let's create a factor with little levels")
#!   x <- ff(letters[4:6], levels=letters[4:6])
#!   message("Let's interpret the same ff file without levels in order to see the codes")
#!   y <- x
#!   levels(y) <- NULL
#!
#!   levels(x)
#!   data.frame(factor=x[], codes=y[])
#!
#!   levels(x) <- appendLevels(levels(x), letters)
#!   levels(x)
#!   data.frame(factor=x[], codes=y[])
#!
#!   x <- sortLevels(x) # implicit recoding is chunked were necessary
#!   levels(x)
#!   data.frame(factor=x[], codes=y[])
#!
#!   message("NEVER forget to reassign the result of recodeLevels or sortLevels, 
#! look at the following mess")
#!   recodeLevels(x, rev(levels(x)))
#!   message("NOW the codings have changed, but not the levels, the result is wrong data")
#!   levels(x)
#!   data.frame(factor=x[], codes=y[])
#!
#!   rm(x);gc()
#!
#! \dontrun{
#!  n <- 5e7
#!
#!  message("reading a factor from a file ist as fast ...")
#!  system.time(
#!  fx <- ff(factor(letters[1:25]), length=n)
#!  )
#!  system.time(x <- fx[])
#!  str(x)
#!  rm(x); gc()
#!
#!
#!  message("... as creating it in-RAM (R-2.11.1) which is theoretically impossible ...")
#!  system.time({
#!  x <- integer(n)
#!  x[] <- 1:25
#!  levels(x) <- letters[1:25]
#!  class(x) <- "factor"
#!  })
#!  str(x)
#!  rm(x); gc()
#!
#!
#!  message("... but is possible if we avoid some  unnecessary copying that is triggered 
#! by assignment functions")
#!  system.time({
#!  x <- integer(n)
#!  x[] <- 1:25
#!  setattr(x, "levels", letters[1:25])
#!  setattr(x, "class", "factor")
#!  })
#!  str(x)
#!  rm(x); gc()
#!
#!  rm(n)
#! }
#!
#! }
#! \keyword{ attribute }


# expand the levels of ... arguments and return joined (unsorted) levels
appendLevels <- function(...){
  unique(unlist(lapply(list(...), function(x){
    if (is.factor(x))
      levels(x)
    else
      x
  })))
}

recodeLevels.factor <- function(x, lev){
  m <- match(levels(x), lev)
  a <- attributes(x)
  a$levels <- lev
  attributes(x) <- NULL
    x <- m[x]
  attributes(x) <- a
  x
}

recodeLevels.ff <- function(x, lev){
  stopifnot(is.factor(x))
  m <- match(levels(x), lev)
  #assigning levels<-NULL kills ramclass
  rc <- ramclass(x)
  levels(x) <- NULL
  i1 <- i2 <- NULL  # keep R CMD CHECK quiet
  ffvecapply({
    i <- hi(i1,i2)
    x[i] <- m[x[i]]
  },X=x)
  #assining levels creates ramclass 'factor' (but not 'ordered')
  levels(x) <- lev
  virtual(x)$ramclass <- rc
  x
}

sortLevels.factor <- function(x){
  l <- levels(x)
  o <- order(l)
  if (identical(o, seq_along(o)))
    x
  else
    recodeLevels.factor(x, l[o])
}

sortLevels.ff <- function(x){
  l <- levels(x)
  o <- order(l)
  if (identical(o, seq_along(o)))
    x
  else
    recodeLevels.ff(x, l[o])
}

sortLevels.ffdf <- function(x){
  for (i in seq_len(ncol(x))){
    y <- x[[i]]
    if (is.factor(y))
      x[[i]] <- sortLevels.ff(y)
  }
  x
}

