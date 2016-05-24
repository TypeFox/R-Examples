# ff coercing of data objects
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-10-09
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/as.ff.R")

#! \name{as.ff}
#! \alias{as.ff}
#! \alias{as.ff.default}
#! \alias{as.ff.ff}
#! \alias{as.ram}
#! \alias{as.ram.default}
#! \alias{as.ram.ff}
#! \title{ Coercing ram to ff and ff to ram objects }
#! \description{
#!    Coercing ram to ff and ff to ram objects while optionally modifying object features.
#! }
#! \usage{
#!   as.ff(x, ...)
#!   as.ram(x, ...)
#!   \method{as.ff}{default}(x, filename = NULL, overwrite = FALSE, ...)
#!   \method{as.ff}{ff}(x, filename = NULL, overwrite = FALSE, ...)
#!   \method{as.ram}{default}(x, ...)
#!   \method{as.ram}{ff}(x, ...)
#! }
#! \arguments{
#!   \item{x}{ any object to be coerced }
#!   \item{filename}{ path and filename }
#!   \item{overwrite}{ TRUE to overwrite the old filename }
#!   \item{\dots}{ \code{\dots} }
#! }
#! \details{
#!   If \command{as.ff.ff} is called on an 'ff' object or \command{as.ram.default} is called on a non-ff object AND no changes are required, the input object 'x' is returned unchanged.
#!   Otherwise the workhorse \code{\link{clone.ff}} is called.
#!   If no change of features are requested, the filename attached to the object remains unchanged, otherwise a new filename is requested (or can be set by the user).
#! }
#! \note{
#!    If you use \code{ram <- as.ram(ff)} for caching, please note that you must \command{\link{close.ff}} before you can write back \code{as.ff(ram, overwrite=TRUE)} (see examples).
#! }
#! \value{
#!   A ram or ff object.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{as.ff.bit}}, \code{\link{ff}}, \code{\link{clone}} %, \code{\link{as.symm}}
#!          , \code{\link{as.vmode}}, \code{\link{vmode}}, \code{\link{as.hi}} }
#! \examples{
#!    message("create ff")
#!    myintff <- ff(1:12)
#!    message("coerce (=clone) integer ff to double ff")
#!    mydoubleff <- as.ff(myintff, vmode="double")
#!    message("cache (=clone) integer ff to integer ram AND close original ff")
#!    myintram <- as.ram(myintff) # filename is retained
#!    close(myintff)
#!    message("modify ram cache and write back (=clone) to ff")
#!    myintram[1] <- -1L
#!    myintff <- as.ff(myintram, overwrite=TRUE)
#!    message("coerce (=clone) integer ram to double ram")
#!    mydoubleram <- as.ram(myintram, vmode="double")
#!    message("coerce (inplace) integer ram to double ram")
#!    myintram <- as.ram(myintram, vmode="double")
#!    message("more classic: coerce (inplace) double ram to integer ram")
#!    vmode(myintram) <- "integer"
#!    rm(myintff, myintram, mydoubleff, mydoubleram); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


as.ram.default <- function(x
, ...  # further parameters to clone.ff (usually none, most notaby 'vmode', if we change something, the filename is lost)
)
{
  if (length(list(...)))
    clone.ff(x              # the workhorse
    , FF_RETURN = FALSE     # tell clone.ff to return a new ram object (old filename is lost)
    , ...
    )
  else
    x
}


as.ram.ff <- function(x
, ... # further paramters to clone.ff (usually none, most notably 'vmode')
){
  if (length(list(...)))
    clone.ff(x                # the workhorse
    , FF_RETURN = FALSE       # tell clone.ff to return a modified ram object (old filename is dropped)
    , ...
    )
  else
    clone.ff(x                # the workhorse
    , FF_RETURN = FALSE       # tell clone.ff to return the identical object as ram
    , filename  = filename(x) # tell clone.ff to keep the filename
    )
}

as.ff.ff <- function(x
, filename  = NULL
, overwrite = FALSE
, ... # further paramters to clone.ff (usually none, most notably 'vmode')
)
{
  if (length(list(...))||!is.null(filename))
    clone.ff(x                # the workhorse
    , FF_RETURN = TRUE        # we require a new ff object, in order to avoid file collisions, we drop the old filename
    , filename  = filename
    , overwrite = overwrite
    , ...
    )
  else
    x
}

as.ff.default <- function(x
, filename  = NULL
, overwrite = FALSE
, ... # further paramters to clone.ff (usually none, most notably 'vmode')
){
  if (length(list(...))||!is.null(filename))
    clone.ff(x                # the workhorse
    , FF_RETURN = TRUE        # we require a new ff object, in order to avoid file collisions, we drop the filename
    , filename  = filename
    , overwrite = overwrite
    , ...
    )
  else
    clone.ff(x                # the workhorse
    , FF_RETURN = TRUE        # tell clone.ff to return the same object as ff
    , filename=filename(x)    # tell clone.ff to keep the name ONLY when no changes
    , overwrite = overwrite
    )
}



if (FALSE){
  library(ff)
  x <- matrix(1:12, 3, 4, dimnames=list(letters[1:3], LETTERS[1:4]))
  y <- as.ff(x)
  z <- as.ram(y)
  a <- as.ff(z)
  close(y);gc()
  a <- as.ff(z)
  a <- as.ff(z, overwrite=TRUE)
  b <- as.ff(z, overwrite=TRUE)


  a <- ff(1:12, dim=c(3,4), dimorder=1:2)
  dimnames(a) <- make.dimnames(a)
  b <- ff(1:12, dim=c(3,4), dimorder=2:1, dimnames=dimnames(a))
  a2 <- ff(1:12, dim=c(3,4), dimorder=1:2, dimnames=dimnames(a), bydim=2:1)
  b2 <- ff(1:12, dim=c(3,4), dimorder=2:1, dimnames=dimnames(a), bydim=2:1)

  names(a) <- 1:length(a)
  names(b) <- 1:length(b)
  names(a2) <- 1:length(a2)
  names(b2) <- 1:length(b2)

  helper <- function(x){
    if (identical(x, a))
      "a"
    else if (identical(x, b))
      "b"
    else if (identical(x, a2))
      "a2"
    else if (identical(x, b2))
      "b2"
  }

  stopifnot(identical(a[], b[]))
  stopifnot(identical(a2[], b2[]))
  stopifnot(identical(a2[], matrix(a[], 3, 4, byrow=TRUE, dimnames=dimnames(a2))))
  stopifnot(identical(b2[], matrix(b[], 3, 4, byrow=TRUE, dimnames=dimnames(b2))))

  for (o in list(a,b,a2,b2))
  for (bydo in list(NULL, 1:2, 2:1))
  for (do in list(NULL, 1:2, 2:1))
  {
    print(list(o=helper(o), do=do, bydo=bydo))
    cc <- clone(o, dimorder=do, bydim=bydo)
    OK <- if (identical(bydo, 2:1)){
      identical(matrix(o[], 3, 4, byrow=TRUE, dimnames=dimnames(o)), cc[]) && identical(names(cc), as.character(as.vector(cc[])))
    }else{
      identical(o[], cc[]) && identical(names(cc), as.character(as.vector(cc[])))
    }
    if (!OK){
      print(cc)
      print(names(cc))
      stop()
    }
  }



}
