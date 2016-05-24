#  R package rjags file R/read.data.R
#  Copyright (C) 2009-2011 Martyn Plummer
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License version
#  2 as published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

read.data <- function(file, format=c("jags","bugs"))
{
    .Deprecated("read.jagsdata", package="rjags")
    format <- match.arg(format)
    switch(format, "jags"=read.jagsdata(file), "bugs"=read.bugsdata(file))
}

read.jagsdata <- function(file)
{
  e <- new.env()
  eval(parse(file), e)
  return(as.list(e))
}

read.bugsdata <- function(file)
{
    bugs.dat <- dget(file)
    for (n in names(bugs.dat)) {
        if (!is.null(dim(bugs.dat[[n]]))) {
            dim(bugs.dat[[n]]) <- rev(dim(bugs.dat[[n]]))
            bugs.dat[[n]] <- aperm(bugs.dat[[n]])
        }
    }
    return(bugs.dat)
}
    

