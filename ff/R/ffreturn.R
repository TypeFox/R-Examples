# Checking ffconform, ffsuitable and creating ffreturn values
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-11-01
# Last changed: 2007-11-14

# source("d:/mwp/eanalysis/ff/R/ffreturn.R")

#! \name{mismatch}
#! \alias{mismatch}
#! \alias{ymismatch}
#! \title{ Test for recycle mismatch }
#! \description{
#!   \command{mismatch} will return TRUE if the larger of nx,ny is not a multiple of the other and the other is >0 (see arithmetic.c).
#!   \command{ymismatch} will return TRUE if nx is not a multiple of ny and ny>0
#! }
#! \usage{
#! mismatch(nx, ny)
#! ymismatch(nx, ny)
#! }
#! \arguments{
#!   \item{nx}{ x length }
#!   \item{ny}{ y length }
#! }
#! \value{
#!   logical scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ffconform}} }
#! \examples{
#!   ymismatch(4,0)
#!   ymismatch(4,2)
#!   ymismatch(4,3)
#!   ymismatch(2,4)
#!   mismatch(4,0)
#!   mismatch(4,2)
#!   mismatch(4,3)
#!   mismatch(2,4)
#! }
#! \keyword{ IO }
#! \keyword{ data }

# symmetric, see arithmetic.c
mismatch <- function(nx,ny){
  if (nx == ny || nx == 1 || ny == 1)
    return(FALSE)
  else if (nx > 0 && ny > 0) {
    if (nx > ny)
      return(as.logical(nx %% ny))
    else
      return(as.logical(ny %% nx))
  }else
    return(FALSE)
}

# asymmetric
ymismatch <- function(nx,ny){
  if (nx == ny || ny == 1)
    return(FALSE)
  else if (ny > 0)
    return(as.logical(nx %% ny))
  else
    return(FALSE)
}


#! \name{dimorderCompatible}
#! \alias{dimorderStandard}
#! \alias{vectorStandard}
#! \alias{dimorderCompatible}
#! \alias{vectorCompatible}
#! \title{ Test for dimorder compatibility }
#! \description{
#!   \command{dimorderStandard} returns TRUE if the dimorder is standard (ascending),
#!   \command{vectorStandard} returns TRUE if the dimorder-bydim combination is compatible with a standard elementwise vector interpretation,
#!   \command{dimorderCompatible} returns TRUE if two dimorders have a compatible  elementwise vector interpretation
#!   and \command{vectorCompatible} returns TRUE if dimorder-bydim combinations have a compatible  elementwise vector interpretation.
#! }
#! \usage{
#! dimorderStandard(dimorder)
#! vectorStandard(dimorder, bydim = NULL)
#! dimorderCompatible(dim, dim2, dimorder, dimorder2)
#! vectorCompatible(dim, dim2, dimorder=NULL, dimorder2=NULL, bydim = NULL, bydim2 = NULL)
#! }
#! \arguments{
#!   \item{dim}{ a \code{\link{dim}} }
#!   \item{dim2}{ a dim }
#!   \item{dimorder}{ a \code{\link{dimorder}} }
#!   \item{dimorder2}{ a dimorder }
#!   \item{bydim}{ a bydim order, see \code{\link{[.ff}} }
#!   \item{bydim2}{ a bydim order, see argument \code{fromdim} in \code{\link{update.ff}} }
#! }
#! \value{
#!   TRUE if compatibility has been detected, FALSE otherwise
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ does not yet gurantee to detect all compatible configurations, but the most important ones }
#! \seealso{ \code{\link{dimorder}}, \code{\link{ffconform}} }
#! \keyword{ IO }
#! \keyword{ data }


# TRUE if the dimorder is standard (or NULL)
dimorderStandard <- function(dimorder)
  is.null(dimorder) || identical(dimorder, 1:length(dimorder))

# TRUE if the vector is standard, FALSE otherwise
vectorStandard <- function(dimorder, bydim=NULL){
  identical(dimorder, bydim) || ( dimorderStandard(dimorder) && dimorderStandard(bydim) )
}

# TRUE if two dim/dimorders are compatible, FALSE otherwise
dimorderCompatible <- function(dim, dim2, dimorder=NULL, dimorder2=NULL){
  if (dimorderStandard(dimorder) && dimorderStandard(dimorder2)){
    return(TRUE)
  }else{
    n <- length(dim)
    n2 <- length(dim2)
    if (n){
      if (n2){
          if (n==n2){
            if (is.null(dimorder))
              dimorder <- 1:n
            if (is.null(dimorder2))
              dimorder2 <- 1:n
            return(identical(dimorder, dimorder2) && identical(dim[dimorder], dim2[dimorder2]))
          }else{
            return(FALSE)
          }
      }else{
        return(dimorderStandard(dimorder))
      }
    }else{
      if (n2){
        return(dimorderStandard(dimorder2))
      }else{
        return(TRUE)
      }
    }
  }
}

# TRUE if dim/dimorder/bydim conflicts with straight vector interpretation
vectorCompatible <- function(
  dim, dim2
, dimorder=NULL, dimorder2=NULL
, bydim=NULL, bydim2=NULL
){
  if(vectorStandard(dimorder, bydim) && vectorStandard(dimorder2, bydim2)){
    return(TRUE)
  }else{
    n <- length(dim)
    n2 <- length(dim2)
    if (n){
      if (n2){
        if (is.null(dimorder))
          dimorder <- 1:n
        if (is.null(dimorder2))
          dimorder2 <- 1:n
        if (is.null(bydim))
          bydim <- 1:n
        if (is.null(bydim2))
          bydim2 <- 1:n
        return(identical(bydim[dimorder], bydim2[dimorder2]) && identical(dim[dimorder], dim2[dimorder2]))
      }else{
        return(vectorStandard(dimorder=dimorder, bydim=bydim))
      }
    }else{
     if (n2){
       return(vectorStandard(dimorder=dimorder2, bydim=bydim2))
     }else{
       return(TRUE)
     }
   }
  }
}


#! \name{ffconform}
#! \alias{ffconform}
#! \title{ Get most conforming argument }
#! \description{
#!   \command{ffconform} returns position of 'most' conformable ff argument or zero if the arguments are not conforming
#! }
#! \usage{
#! ffconform(\dots, vmode = NULL, fail = "stop")
#! }
#! \arguments{
#!   \item{\dots}{ two or more ff objects }
#!   \item{vmode}{ handing over target vmode here supresses searching for a common vmode, see \code{\link{maxffmode}} }
#!   \item{fail}{ the name of a function to call if not-conforming, default \code{\link{stop}} }
#! }
#! \details{
#!   A reference argument is defined to be the first argument with a \code{\link[ff:dim.ff]{dim}} attribute or the longest vector.
#!   The other arguements are then compared to the reference to check for conformity,
#!   which is violated if vmodes are not conforming
#!   or if the reference has not a multiple length of each other
#!   or if the dimensions do not match
#!   or if we have a dimorder conflict because not all arguments have the same \code{\link{dimorderStandard}}.
#! }
#! \value{
#!   the position of the most conforming argument or 0 (zero) if not conforming.
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ xx Work in progress for package \pkg{R.ff} }
#! \seealso{ \code{\link{ffsuitable}}, \code{\link{maxffmode}}, \code{\link{ymismatch}}, \code{\link{stop}}, \code{\link{warning}}, \code{\link{dimorderStandard}} }
#! \examples{
#!   a <- ff(1:10)
#!   b <- clone(a)
#!   c <- ff(1:20)
#!   d <- ff(1:21)
#!   ffconform(a,b)
#!   ffconform(c,a)
#!   ffconform(a,c)
#!   ffconform(c,a,b)
#!
#!   d1 <- ff(1:20, dim=c(2,10))
#!   d2 <- ff(1:20, dim=c(10,2))
#!   ffconform(c,d1)
#!   ffconform(c,d2)
#!   ffconform(d1,c)
#!   ffconform(d2,c)
#!   try(ffconform(d1,d2))
#!   ffconform(d1,d1)
#!
#!   rm(a,b,c,d1,d2); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

ffconform <- function(
  ...
, vmode = NULL
, fail = "stop"
)
{
  l <- list(...)
  nl <- length(l)

  # determine smallest lossless vmode if not given
  if (is.null(vmode)){
    vmode <- lapply(l,vmode)
    ffmode <- do.call("maxffmode", vmode)
    if (ffmode){
      vmode <- names(ffmode)
    }else{
      if (!is.null(fail)) do.call(fail, list("no common vmode"))
      return(0L)
    }
  }

  # determine reference: first object with dim or largest vector
  a <- 1L
  if (nl>1){
    ff <- l[[1L]]
    if (is.null(dim(ff))){
      n <- length(ff)
      for (i in 2:nl){
        ff2 <- l[[i]]
        if (!is.null(dim(ff2))){
          a <- i
          break
        }
        n2 <- length(ff2)
        if (n2>n){
          a <- i
          n <- n2
        }
      }
    }
    ff <- l[[a]]
    n <- length(ff)
    d <- dim(ff)
    do <- dimorder(ff)
    don <- !dimorderStandard(do)
    sym <- symmetric(ff)
    dia <- fixdiag(ff)
    # check conformity with reference
    for (i in (1:nl)[-a]){
      ff2 <- l[[i]]
      n2 <- length(ff2)
      d2 <- dim(ff2)
      if (is.null(d2)){
        if (ymismatch(n,n2)) {
          if (!is.null(fail)) do.call(fail, list("longer object is not multiple of shorter one"))
          return(0L)
        }
        if (don) {
          if (!is.null(fail)) do.call(fail, list("trying to mix vectors and array with non-standard dimorder"))
          return(0L)
        }
      }else{
        if (!identical(d,d2)){
          if (!is.null(fail)) do.call(fail, list("non-conforming arrays"))
          return(0L)
        }
        do2 <- dimorder(ff2)
        if ( !identical(do,do2) && ( don || !dimorderStandard(do2) )){
          if (!is.null(fail)) do.call(fail, list("arrays with conflicting dimorder"))
          return(0L)
        }
        if (!identical(sym, symmetric(ff2))){
          if (!is.null(fail)) do.call(fail, list("mixing symmetric and non-symmetric"))
          return(0L)
        }
        if (!identical(dia, fixdiag(ff2))){
          if (!is.null(fail)) do.call(fail, list("mixing fixed and non-fixed diagonal"))
          return(0L)
        }
      }
    }
  }
  attr(a, "vmode") <- vmode
  return(a)
}


#! \name{ffsuitable}
#! \alias{ffsuitable}
#! \alias{ffsuitable_attribs}
#! \title{ Test ff object for suitability }
#! \description{
#!   \command{ffsuitable} tests whether \code{FF_RETURN} is an \code{\link{ff}} object like \code{FF_PROTO} and having attributes \code{FF_ATTR}.
#! }
#! \usage{
#! ffsuitable(FF_RETURN, FF_PROTO = NULL, FF_ATTR = list()
#! , strict.dimorder = TRUE, fail = "warning")
#! ffsuitable_attribs(x)
#! }
#! \arguments{
#!   \item{x}{ an object from which to extract attributes for comparison }
#!   \item{FF_RETURN}{ the object to be tested for suitability }
#!   \item{FF_PROTO}{ the prototype object which \code{FF_RETURN} should match }
#!   \item{FF_ATTR}{ a list of additional attributes dominating those from \code{FF_PROTO} }
#!   \item{strict.dimorder}{ if TRUE ffsuitability requires that the dimorders are standard (ascending) }
#!   \item{fail}{ name of a function to be called if not ffsuitable (default \code{\link{warning}}) }
#! }
#! \value{
#!   TRUE if \code{FF_RETURN} object is suitable, FALSE otherwise
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ xx Work in progress for package \pkg{R.ff} }
#! \seealso{ \code{\link{ffconform}}, \code{\link{ffreturn}} }
#! \keyword{ IO }
#! \keyword{ data }


ffsuitable_attribs <- function(x){
  list(
    vmode       = vmode(x)
  , length      = length(x)
  , dim         = dim(x)
  , dimorder    = dimorder(x)
  , symmetric   = symmetric(x)
  , fixdiag     = fixdiag(x)
  )
}

# check whether x is suitable to absorb an ff result with given specification
# returns TRUE or FALSE and optionally warns
ffsuitable <- function(
   FF_RETURN
 , FF_PROTO = NULL          # from FF_PROTO compares vmode, length, dim, dimorder, symmetric, fixdiag
 , FF_ATTR = list()  # compares list of explicit attributes for comparing. NOTE the spelling of the attributes above ('dim', not 'Dim'). NOTE that length identical <number>L, not identical <number>
 , strict.dimorder = TRUE
 , fail = "warning"
)
{
  if (!inherits(FF_RETURN, "ff"))
    return(FALSE)
  xlist <- ffsuitable_attribs(FF_RETURN)
  if (is.null(FF_PROTO)){
    ylist <- FF_ATTR
  }else{
    ylist <- ffsuitable_attribs(FF_PROTO)
    if (length(FF_ATTR))
      ylist[names(FF_ATTR)] <- FF_ATTR
  }

  # first check vector conditions
  if (is.null(ylist$vmode) || ylist$vmode=="NULL")
    stop("need vmode")
  if (!identical(as.vector(xlist$vmode), as.vector(ylist$vmode))){
    if (!is.null(fail))
      do.call(fail, list(paste(deparse(substitute), "not suitable because vmode differs")))
    return(FALSE)
  }
  if (is.null(ylist$length)){
    if (is.null(ylist$dim))
      stop("need length or dim")
    else
      ylist$length <- as.integer(prod(ylist$dim))
  }
  if (!identical(as.vector(xlist$length), as.vector(ylist$length))){
    if (!is.null(fail))
      do.call(fail, list(paste(deparse(substitute), "not suitable because length differs")))
    return(FALSE)
  }
  if (is.null(ylist$dim)){
    if (is.null(xlist$dim)){
      return(TRUE)
    }else{
      if (strict.dimorder && !identical(xlist$dimorder, sort(xlist$dimorder))){
        if (!is.null(fail))
          do.call(fail, list(paste(deparse(substitute), "not suitable because is has non-standard dimorder")))
        return(FALSE)
      }
    }
  }else{
    if (is.null(xlist$dim)){
      if (strict.dimorder){
        if (is.null(ylist$dimorder))
          stop("dimorder not provided for strict.dimorder")
        if (!identical(ylist$dimorder, sort(ylist$dimorder))){
          if (!is.null(fail))
            do.call(fail, list(paste(deparse(substitute), "not suitable because strict non-standard dimorder was required")))
          return(FALSE)
        }
      }
      return(TRUE)
    }else{
      if (!identical(as.vector(xlist$dim), as.vector(ylist$dim))){
        if (!is.null(fail))
          do.call(fail, list(paste(deparse(substitute), "not suitable because of differing dimensions")))
        return(FALSE)
      }
      if (!length(ylist[["symmetric"]])){
        ylist$symmetric <- FALSE
      }
      if (!length(ylist[["fixdiag"]])){
        ylist <- c(ylist, list(fixdiag=NULL))
      }
      if (!identical(xlist[c("symmetric","fixdiag")], ylist[c("symmetric","fixdiag")])){
        if (!is.null(fail))
          do.call(fail, list(paste(deparse(substitute), "not suitable because of differing matrix types (w/r to symmetric, fixdiag)")))
        return(FALSE)
      }
    }
  }
  TRUE
}



#! \name{ffreturn}
#! \alias{ffreturn}
#! \title{ Return suitable ff object }
#! \description{
#!   \command{ffreturn} returns \code{FF_RETURN} if it is \code{\link{ffsuitable}} otherwise creates a suitable \code{\link{ffsuitable}} object
#! }
#! \usage{
#! ffreturn(FF_RETURN = NULL, FF_PROTO = NULL, FF_ATTR = NULL)
#! }
#! \arguments{
#!   \item{FF_RETURN}{ the object to be tested for suitability }
#!   \item{FF_PROTO}{ the prototype object which \code{FF_RETURN} should match }
#!   \item{FF_ATTR}{ a list of additional attributes dominating those from \code{FF_PROTO} }
#! }
#! \value{
#!   a suitable \code{\link{ffsuitable}} object
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ xx Work in progress for package \pkg{R.ff} }
#! \seealso{ \code{\link{ffconform}}, \code{\link{ffsuitable}} }
#! \keyword{ IO }
#! \keyword{ data }


ffreturn <- function(FF_RETURN=NULL, FF_PROTO=NULL, FF_ATTR=NULL){
  if (is.null(FF_RETURN)){
    if (is.null(FF_PROTO))
      FF_RETURN <- TRUE
    else
      FF_RETURN <- is.ff(FF_PROTO)
  }
  if (ffsuitable(FF_RETURN, FF_PROTO, FF_ATTR=FF_ATTR)){
    if (!is.null(FF_ATTR[["initdata"]]))
      FF_RETURN[] <- FF_ATTR[["initdata"]]
    FF_RETURN
  }else{
    if (!is.logical(FF_RETURN))
      FF_RETURN <- TRUE
    if (is.null(FF_PROTO))
      do.call("ff", c(FF_ATTR, list(FF_RETURN=FF_RETURN)))  # ff will create whatever is needed
    else
      do.call("clone.ff", c(list(FF_PROTO), FF_ATTR, list(FF_RETURN=FF_RETURN)))  # clone.ff will create whatever is needed
  }
}
