# clone utilities for bit,bit64,ff
# (c) 2014 Jens Oehlschlägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2014-03-02



#! \name{clone}
#! \alias{clone}
#! \alias{clone.list}
#! \alias{clone.default}
#! \alias{still.identical}
#! \title{ Cloning ff and ram objects }
#! \description{
#!   \command{clone} physically duplicates objects and can additionally change some features, e.g. length.
#! }
#! \usage{
#! clone(x, \dots)
#! \method{clone}{list}(x, \dots)
#! \method{clone}{default}(x, \dots)
#! still.identical(x, y)
#! }
#! \arguments{
#!   \item{x}{ \code{x} }
#!   \item{y}{ \code{y} }
#!   \item{\dots}{ further arguments to the generic }
#! }
#! \details{
#!   \command{clone} is generic. 
#!   \command{clone.default} currently only handles atomics. 
#!   \command{clone.list} recursively clones list elements.
#!   \command{still.identical} returns TRUE if the two atomic arguments still point to the same memory.
#! }
#! \value{
#!   an object that is a deep copy of x
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[ff]{clone.ff}} }
#! \examples{
#!   x <- 1:12
#!   y <- x
#!   still.identical(x,y)
#!   y[1] <- y[1]
#!   still.identical(x,y)
#!   y <- clone(x)
#!   still.identical(x,y)
#!   rm(x,y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

# named <- function(x)
  # .Call("R_bit_named", x, PACKAGE="bit")

still.identical <- function(x, y){
  .Call("r_ram_truly_identical", x = x, y = y, PACKAGE = "bit")
}

clone.default <- function(x
, ... # passed to clone
){
  if (is.atomic(x)){
    if (length(x))
      x[1] <- x[1]  # force a copy around COPY ON MODIFY
    x
  }else{
    stop("clone not defined for type")
  }
}

clone.list <- function(x
, ... # passed to clone
){
  lapply(x, clone, ...)
}

