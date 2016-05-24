# Utilities for ff
# (c) 2007 Jens Oehlschlägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-08-24
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/util.R")


#! \name{unclass_-}
#! \alias{unclass<-}
#! \title{ Unclassed assignement }
#! \description{
#!   With \command{unclass<-} you can circumvent class dispatch on the assignment operator
#! }
#! \usage{
#! unclass(x) <- value
#! }
#! \arguments{
#!   \item{x}{ some object }
#!   \item{value}{ the value to be assigned }
#! }
#! \value{
#!   the modified object
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{unclass}}, \code{\link{undim}} }
#! \examples{
#!   x <- factor(letters)
#!   unclass(x)[1:3] <- 1L
#!   x
#! }
#! \keyword{ IO }
#! \keyword{ data }


"unclass<-" <- function(x, value){
  cl <- class(x)
  x <- unclass(x)
  x <- value
  class(x) <- cl
  x
}


#! \name{undim}
#! \alias{undim}
#! \title{ Undim }
#! \description{
#!   \command{undim} returns its input with the dim attribute removed
#! }
#! \usage{
#! undim(x)
#! }
#! \arguments{
#!   \item{x}{ an object }
#! }
#! \value{
#!   x without dim attribute
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{unclass<-}}, \code{\link{unclass}}, \code{\link{unname}}, \code{\link{dim}} }
#! \examples{
#!   x <- matrix(1:12, 3)
#!   x
#!   undim(x)
#! }
#! \keyword{ IO }
#! \keyword{ data }

undim <- function(x){
  dim(x) <- NULL
  x
}



#! \name{repnam}
#! \Rdversion{1.1}
#! \alias{repnam}
#! \title{
#!   Replicate with names
#! }
#! \description{
#!   Function \code{repnam} replicates its \code{argument} to the desired \code{length}, either by simply \code{\link{rep}licating} or - if it has \code{\link{names}} - by replicating the \code{default} and matching the argument by its names.
#! }
#! \usage{
#! repnam(argument, names = NULL, len=length(names), default = list(NULL))
#! }
#! \arguments{
#!   \item{argument}{
#!   a named or non-named vector or list to be replicated
#! }
#!   \item{names}{
#!   NULL or a charcter vector of names to which the argument names are matched
#! }
#!   \item{len}{
#!   the desired length (required if names is not given)
#! }
#!   \item{default}{
#!   the desired default which is replicated in case names are used (the default \code{list(NULL)} is suitable for a list argument)
#! }
#! }
#! \value{ an object like argument or default having length len }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \note{
#!   This is for internal use, e.g. to handle argument \code{colClasses} in \code{\link{read.table.ffdf}}
#! }
#! \seealso{
#!   \code{\link{rep}}, \code{\link{vector}}, \code{\link[bit]{repfromto}}
#! }
#! \examples{
#!  message("a list example")
#!  repnam(list(y=c(1,2), z=3), letters)
#!  repnam(list(c(1,2), 3), letters)
#!
#!  message("a vector example")
#!  repnam(c(y=1, z=3), letters, default=NA)
#!  repnam(c(1, 3), letters, default=NA)
#!
#! }
#! \keyword{ utilities }

repnam <- function(argument, names=NULL, len=length(names), default=list(NULL)){
  argnam <- names(argument)
  if (is.null(argnam)){
    ret <- rep(argument, length.out=len)
  }else{
    if (is.null(names))
      stop("cannot match named argument to non-named prototype")
    i <- match(argnam, names)
    if (any(is.na(i)))
      stop("the following argument names do not match", paste("'", argnam[is.na(i)], "'", sep="", collapse=","))
    ret <- rep(default, len)
    names(ret) <- names
    ret[i] <- argument
  }
  ret
}
