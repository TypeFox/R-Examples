# Generic functions for ff
# (c) 2007 Daniel Adler
# Licence: free BSD
# Provided 'as is', use at your own risk
# Created: 2007-08-24
# Last changed: 2008-07-26

#! \name{getpagesize}
#! \alias{getpagesize}
#! \alias{getdefaultpagesize}
#! \alias{getalignedpagesize}
#! \title{
#!   Get page size information
#!   }
#! \description{
#!   The function is used for obtaining the natural OS-specific page size in Bytes.
#!   \command{getpagesize} returns the OS-specific page size in Bytes for memory mapped files, while \command{getdefaultpagesize} returns a suggested page size.
#!   \command{getalignedpagesize} returns the pagesize as a multiple of the OS-specific page size in Bytes, which is the correct way to specify pagesize in ff.
#!   }
#! \usage{
#!   getpagesize()
#!   getdefaultpagesize()
#!   getalignedpagesize(pagesize)
#!  }
#! \arguments{
#!   \item{pagesize}{ a desired pagesize in bytes }
#!   }
#! \value{
#!   An integer giving the page size in Bytes.
#! }
#! \author{ Daniel Adler, Jens Oehlschlägel }
#! \examples{
#!   getpagesize()
#!   getdefaultpagesize()
#!   getalignedpagesize(2000000)
#!   }
#! \keyword{IO}

getpagesize <- function()
  .Call("getpagesize", PACKAGE="ff")

getdefaultpagesize <- function()
  #if (.Platform$OS=="windows")
    getalignedpagesize(65536L)
  #else
  #  getalignedpagesize(1048576L)

getalignedpagesize <- function(pagesize)
  if (pagesize){
    syspagesize <- getpagesize()
    if (pagesize %% syspagesize){
      alignedpagesize <- as.integer((pagesize %/% syspagesize + 1L) * syspagesize)
      warning("aligning pagesize ", pagesize, " with systempagesize ", syspagesize, " to ", alignedpagesize)
      alignedpagesize
    }else
      as.integer(pagesize)
  }else{
    stop("rethink this")
  }

