# Matrix transpose and virtual matrix transpose for ff objects
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-11-28
# Last changed: 2007-11-28

# source("d:/mwp/eanalysis/ff/R/vt.R")


#! \name{vt}
#! \alias{vt}
#! \alias{vt.ff}
#! \alias{vt.default}
#! \alias{t.ff}
#! \title{ Virtual transpose }
#! \description{
#!   The \command{vt} generic does a matrix  or array transpose by modifying \code{\link[=virtual.ff]{virtual}} attributes
#!   rather than by physically copying matrix elements.
#! }
#! \usage{
#! vt(x, \dots)
#! \method{vt}{ff}(x, \dots)
#! \method{vt}{default}(x, \dots)
#! \method{t}{ff}(x)
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   The \code{vt.ff} method does transpose through reversing \code{\link{dim.ff}} and \code{\link{dimorder}}.
#!   The \code{vt.default} method is a wrapper to the standard transpose \code{\link{t}}. \cr
#!   The \code{t.ff} method creates a transposed \code{\link{clone}}. \cr
#!   If \code{x} has a virtual window \code{\link{vw}} defined, \code{vt.ff} returns an ff object with a transposed virtual window,
#!   the \code{t.ff} method return a transposed clone of the virtual window content only.
#! }
#! \value{
#!   an object that behaves like a transposed matrix
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{dim.ff}}, \code{\link{vw}}, \code{\link[=virtual.ff]{virtual}} }
#! \examples{
#!   x <- ff(1:20, dim=c(4,5))
#!   x
#!   vt(x)
#!   y <- t(x)
#!   y
#!   vw(x) <- cbind(c(1,3,0),c(1,4,0))
#!   x
#!   vt(x)
#!   y <- t(x)
#!   y
#!   rm(x,y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ array }

vt.default <- function(x
, ... # dummy to keep R CMD check quiet
)
  t(x)

vt.ff <- function(x
, ... # dummy to keep R CMD check quiet
){
  if (symmetric(x))
    return(x)
  d <- dim(x)
  if (is.null(d) )
    stop("not an array")
  if (length(d)!=2)
    warning("not a matrix")
  vw <- vw(x)
  if (is.null(vw)){
    do <- dimorder(x)
    dn <- dimnames(x)
    dim(x) <- rev(d)
    dimorder(x) <- rev(do)
    dimnames(x) <- rev(dn)
  }else{
    vw(x) <- NULL
    d <- dim(x)
    do <- dimorder(x)
    dn <- dimnames(x)
    dim(x) <- d[2:1]
    dimorder(x) <- rev(do)
    dimnames(x) <- rev(dn)
    vw(x) <- vw[,ncol(vw):1,drop=FALSE]
  }
  x
}

t.ff <- function(x){
  if (symmetric(x))
    return(x)
  clone(vt(x), dimorder=rev(dimorder(x)))
}
