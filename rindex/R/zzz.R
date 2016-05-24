.onLoad <- function(lib, pkg) {
       warning("package rindex is under development and still may change considerably")
       library.dynam("rindex", pkg, lib)
}

.onUnload <- function(libpath){
       library.dynam.unload("rindex", libpath)
}


#! \name{indexInit}
#! \alias{indexInit}
#! \alias{indexDone}
#! \title{ Load / unload rindex library }
#! \description{
#!   Function \code{indexInit} loads the rindex shared library, \code{indexDone} unloads the rindex shared library.
#! }
#! \usage{
#! indexInit()
#! indexDone()
#! }
#! \details{
#!   You are responsible to free all memory using \code{\link{indexDelTree}} before calling \code{indexDone}.
#! }
#! \value{
#!   See \code{\link[base]{dyn.load}} and \code{\link[base]{dyn.unload}}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{dyn.load}}, \code{\link[base]{dyn.unload}}, \code{\link{indexDelTree}} }
#! \keyword{ misc }
#! \keyword{ database }


# library level
indexInit <- function(){
  dyn.load(file.path(.libPaths(), "rindex", "libs", paste("rindex", .Platform$dynlib.ext, sep = "")))
}
indexDone <- function(){
  dyn.unload(file.path(.libPaths(), "rindex", "libs", paste("rindex", .Platform$dynlib.ext, sep = "")))
}


#! \name{indexDemoClose}
#! \alias{indexDemoOpen}
#! \alias{indexDemoClose}
#! \title{ Demo functions for creating and removing external pointers }
#! \description{
#!   Function \code{indexDemoOpen} creates an external pointer pointing to an integer vector of length \option{n}, \code{indexDemoClose} frees the memory linked to the pointer.
#! }
#! \usage{
#! indexDemoOpen(n = 10)
#! indexDemoClose(extPtr)
#! }
#! \arguments{
#!   \item{extPtr}{ the external pointer returned by \code{indexDemoOpen} }
#!   \item{n}{ the length of the vector stored in C RAM }
#! }
#! \details{
#!   If the returned pointer is removed finalization happens at the next \code{gc()}, if \code{indexDemoClose} is used on the pointer then finalization happens immediately, if the pointer is then removed, the finalizer is called \emph{again} at the next \code{gc()} (but doesn't finalize again).
#! }
#! \value{
#!   Function \code{indexDemoOpen} returns an external pointer, \code{indexDemoClose} returns the integer vector previously linked by the pointer.
#! }
#! \references{ R Development Core Team (2007). Writing R Extensions. }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{indexInit}}, \code{\link{indexDone}}, \code{\link{indexAddTree}}, \code{\link{indexDelTree}} }
#! \examples{
#! ptr <- indexDemoOpen()
#! rm(ptr)
#! gc()
#!
#! ptr <- indexDemoOpen()
#! indexDemoClose(ptr)
#! rm(ptr)
#! gc()
#! }
#! \keyword{ misc }


# test external pointers
indexDemoOpen <- function(n=10){
  .Call("demo_rindex_open", as.integer(n), PACKAGE="rindex")
}
indexDemoClose <- function(extPtr){
  .Call("demo_rindex_close", extPtr, PACKAGE="rindex")
}


#! \name{strncmp}
#! \alias{strcmp}
#! \alias{strncmp}
#! \title{ String comparison á la C }
#! \description{
#!   Functiions to compare two vectors of strings like the standard C functions do.
#! }
#! \usage{
#! strcmp(a, b)
#! strncmp(a, b, n)
#! }
#! \arguments{
#!   \item{a}{ character vector }
#!   \item{b}{ character vector }
#!   \item{n}{ number of characters to compare }
#! }
#! \value{
#!   1 if a>b, -1 if a<b, 0 if a==b.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{substr}}, \code{\link[base]{Comparison}} }
#! \keyword{ misc }



# emulating C strcmp but not used
strcmp <- function(a,b){
  ifelse(a>b, 1, ifelse(a<b, -1, 0))
}
strncmp <- function(a,b,n){
  a <- substr(a,1,n)
  b <- substr(b,1,n)
  ifelse(a>b, 1, ifelse(a<b, -1, 0))
}
