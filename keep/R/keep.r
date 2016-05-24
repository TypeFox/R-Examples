#  Part of the R package keep,
#  Arrays with better control over dimension dropping
#  --------------------------------------------------
#
#  By Paavo Jumppanen
#  Copyright (C) 2015 CSIRO
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


pkg.env <- new.env()
pkg.env$keep.flag <- FALSE


keep<-function(index)
{
  pkg.env$keep.flag = TRUE
  return (index)
}


kOarray <- function(data=NA, dim=length(data), dimnames=NULL, offset=rep(1, length(dim)), drop.negative=TRUE)
{
  if (requireNamespace("Oarray", quietly = TRUE))
  {
    res <- Oarray::Oarray(data, dim, dimnames, offset, drop.negative)
    class(res) <- c("keep", class(res))
    return(res)
  }
  else
  {
    stop("Cannot load package Oarray")
  }
}


karray <- function(data = NA, dim = length(data), dimnames = NULL)
{
  res <- array(data, dim, dimnames)
  class(res) <- c("keep", class(res))
  return(res)
}


"as.karray" <- function(x, ...)
{
  x <- as.array(x)
  class(x) <- c("keep", class(x))
  return(x)
}


"as.kOarray" <- function(x, offset=rep(1, length(dim)), drop.negative=TRUE)
{
  if (requireNamespace("Oarray", quietly = TRUE))
  {
    x <- Oarray::as.Oarray(x)
    class(x) <- c("keep", class(x))
    return(x)
  }
  else
  {
    stop("Cannot load package Oarray")
  }
}


"as.array.keep" <- function(x, ...)
{
  x <- unclass(x)
  NextMethod(x)
}


"[.keep" <- function(x, ...)
{
  arguments <- as.list(substitute(list(...)))[-1L]
  src.dims <- attr(x,"dim")

  if (length(arguments) > length(src.dims))
  {
    stop('too many arguments')
  }
  else if (length(arguments) < length(src.dims))
  {
    # Too few arguments but assume indexing arrays and let the
    # default processing handle it.
    res <- NextMethod()
  }
  else
  {
    dest.dims <- c()

    for (cn in 1:length(arguments))
    {
      arg.expr <- deparse(arguments[[cn]])
      null.arg <- (arg.expr == '')

      if (null.arg)
      {
        dest.dims <- append(dest.dims, src.dims[cn])
      }
      else
      {
        pkg.env$keep.flag <- FALSE
        indices <- eval(arguments[[cn]], envir=parent.frame())
        keep.flag <- pkg.env$keep.flag
        dimension <- length(indices)

        if ((dimension > 1) || keep.flag)
          dest.dims <- append(dest.dims, dimension)
      }
    }

    res <- NextMethod()

    if (length(dest.dims) > 1)
    {
      if (class(res)[1] != "keep")
        class(res) <- c("keep", class(res))

      dim(res) <- dest.dims
    }
  }

  return (res)
}


# Add S4 class unions so we can use these arrays in S4 classes
setOldClass(c("keep", "Oarray"))
setClassUnion("karray", c("keep", "array"))
setClassUnion("kOarray", c("keep", "Oarray"))
