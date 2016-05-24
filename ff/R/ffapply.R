# Batch utilities for operations on ff objects
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-09-03
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/ffapply.R")



#! \name{ffapply}
#! \alias{ffapply}
#! \alias{ffvecapply}
#! \alias{ffrowapply}
#! \alias{ffcolapply}
#! \title{ Apply for ff objects }
#! \description{
#!   The \code{ffapply} functions support convenient batched processing of ff objects
#!   such that each single batch or chunk will not exhaust RAM
#!   and such that batchs have sizes as similar as possible, see \code{\link[bit]{bbatch}}.
#!   Differing from R's standard \code{\link{apply}} which applies a \code{FUNction},
#!   the \code{ffapply} functions do apply an \code{EXPRession} and provide two indices \code{FROM="i1"} and \code{TO="i2"},
#!   which mark beginning and end of the batch and can be used in the applied expression.
#!   Note that the ffapply functions change the two indices in their parent frame, to avoid conflicts you can use different names through \code{FROM="i1"} and \code{TO="i2"}.
#!   For support of creating return values see details.
#! }
#! \usage{
#! ffvecapply(EXPR, X = NULL, N = NULL, VMODE = NULL, VBYTES = NULL, RETURN = FALSE
#! , CFUN = NULL, USE.NAMES = TRUE, FF_RETURN = TRUE, BREAK = ".break"
#! , FROM = "i1", TO = "i2"
#! , BATCHSIZE = .Machine$integer.max, BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE)
#! ffrowapply(EXPR, X = NULL, N = NULL, NCOL = NULL, VMODE = NULL, VBYTES = NULL
#! , RETURN = FALSE, RETCOL = NCOL, CFUN = NULL, USE.NAMES = TRUE, FF_RETURN = TRUE
#! , FROM = "i1", TO = "i2"
#! , BATCHSIZE = .Machine$integer.max, BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE)
#! ffcolapply(EXPR, X = NULL, N = NULL, NROW = NULL, VMODE = NULL, VBYTES = NULL
#! , RETURN = FALSE, RETROW = NROW, CFUN = NULL, USE.NAMES = TRUE, FF_RETURN = TRUE
#! , FROM = "i1", TO = "i2"
#! , BATCHSIZE = .Machine$integer.max, BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE)
#! ffapply(EXPR = NULL, AFUN = NULL, MARGIN = NULL, X = NULL, N = NULL, DIM = NULL
#! , VMODE = NULL, VBYTES = NULL, RETURN = FALSE, CFUN = NULL, USE.NAMES = TRUE
#! , FF_RETURN = TRUE, IDIM = "idim"
#! , FROM = "i1", TO = "i2", BREAK = ".break"
#! , BATCHSIZE = .Machine$integer.max, BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{EXPR}{ the \code{\link{expression}} to be applied }
#!   \item{AFUN}{ \code{ffapply} only: alternatively to \code{EXPR} the name of a function to be applied, automatically converted to \code{EXPR} }
#!   \item{MARGIN}{ \code{ffapply} only: the margins along which to loop in \code{ffapply} }
#!
#!   \item{X}{ an ff object from which several parameters can be derived, if they are not given directly: \code{N, NCOL, NROW, DIM, VMODE, VBYTES, FF_RETURN} }
#!
#!   \item{N}{ the total number of elements in the loop, e.g. number of elements in \code{ffvecapply} or number of rows in \code{ffrowapply} }
#!   \item{NCOL}{ \code{ffrowapply} only: the number of columns needed to calculate batch sizes }
#!   \item{NROW}{ \code{ffcolapply} only: the number of rows needed to calculate batch sizes }
#!   \item{DIM}{ \code{ffapply} only: the dimension of the array needed to calculate batch sizes }
#!
#!   \item{VMODE}{ the \code{\link{vmode}} needed to prepare the \code{RETURN} object and to derive \code{VBYTES} if they are not given directly }
#!   \item{VBYTES}{ the bytes per cell -- see \code{\link{.rambytes}} -- to calculate the RAM requirements per cell }
#!   \item{BATCHBYTES}{ the max number of bytes per batch, default \code{getOption("ffbatchbytes")} }
#!   \item{BATCHSIZE}{ an additional restriction on the number of loop elements, default=\code{.Machine$integer.max} }
#!
#!   \item{FROM}{ the name of the index that marks the beginning of the batch, default 'i1', change if needed to avoid naming-conflicts in the calling frame }
#!   \item{TO}{ the name of the index that marks the end of the batch, default 'i2', change if needed to avoid naming-conflicts in the calling frame }
#!   \item{IDIM}{ \code{ffapply} only: the name of an R variable used for loop-switching, change if needed to avoid naming-conflicts in the calling frame }
#!   \item{BREAK}{ \code{ffapply} only: the name of an R object in the calling frame that triggers break out of the batch loop, if 1) it exists 2) is.logical and 3) is TRUE }
#!
#!   \item{RETURN}{ \code{TRUE} to prepare a return value (default \code{FALSE}) }
#!   \item{CFUN}{ name of a collapsing function, see \code{\link{CFUN}} }
#!   \item{RETCOL}{ \code{NULL} gives return \code{vector[1:N]}, \code{RETCOL} gives return \code{matrix[1:N, 1:RETCOL]} }
#!   \item{RETROW}{ \code{NULL} gives return \code{vector[1:N]}, \code{RETROW} gives return \code{matrix[1:RETROW, 1:N]} }
#!   \item{FF_RETURN}{ \code{FALSE} to return a ram object, \code{TRUE} to return an ff object, or an ff object that is \code{\link{ffsuitable}} to absorb the return data }
#!   \item{USE.NAMES}{ \code{FALSE} to suppress attaching names or dimnames to the result }
#!
#!   \item{VERBOSE}{ \code{TRUE} to verbose the batches }
#! }
#! \details{
#!   \command{ffvecapply} is the simplest ffapply method for \code{ff_vectors}. \command{ffrowapply} and \command{ffcolapply} is for \code{ff_matrix},
#!   and \command{ffapply} is the most general method for \code{ff_array}s and \code{ff_vector}s.
#!   \cr
#!   There are many ways to change the return value of the ffapply functions.
#!   In its simplest usage -- batched looping over an expression -- they don't return anything, see \code{\link{invisible}}.
#!   If you switch \code{RETURN=TRUE} in \command{ffvecapply} then it is assumed that all looped expressions together return one vector of length \code{N},
#!   and via parameter \code{FF_RETURN}, you can decide whether this vector is in ram or is an ff object (or even which ff object to use).
#!   \command{ffrowapply} and \command{ffcolapply} additionally have parameter \code{RETCOL} resp. \code{RETROW} which defaults to returning a matrix of the original size;
#!   in order to just return a vector of length \code{N} set this to \code{NULL}, or specify a number of columns/rows for the return matrix.
#!   It is assumed that the expression will return appropriate pieces for this return structure (see examples).
#!   If you specify \code{RETURN=TRUE} and a collapsing function name \code{CFUN}, then it is assumed that the batched expressions return aggregated information,
#!   which is first collected in a list, and finally the collapsing function is called on this list: \code{do.call(CFUN, list)}. If you want to return the unmodified list,
#!   you have to specify \code{CFUN="list"} for obvious reasons.
#!   \cr
#!   \code{ffapply} allows usages not completly unlike \code{\link{apply}}: you can specify the name of a function \code{AFUN} to be applied over \code{MARGIN}.
#!   However note that you must specify \code{RETURN=TRUE} in order to get a return value.
#!   Also note that currently ffapply assumes that your expression returns exactly one value per cell in \code{DIM[MARGINS]}.
#!   If you want to return something more complicated, you MUST specify a \code{CFUN="list"} and your return value will be a list with dim attribute \code{DIM[MARGINS]}.
#!   This means that for a function \code{AFUN} returning a scalar, \code{ffapply} behaves very similar to \code{\link{apply}}, see examples.
#!   Note also that \code{ffapply} might create a object named '.ffapply.dimexhausted' in its parent frame,
#!   and it uses a variable in the parent frame for loop-switching between dimensions, the default name 'idim' can be changed using the \code{IDIM} parameter.
#!   Finally you can break out of the implied loops by assigning \code{TRUE} to a variable with the name in \code{BREAK}.
#! }
#! \value{
#!   see details
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ xx The complete generation of the return value is preliminary and the arguments related to defining the return value might still change, especially ffapply is work in progress }
#! \seealso{ \code{\link{apply}}, \code{\link{expression}}, \code{\link[bit]{bbatch}}, \code{\link[bit]{repfromto}}, \code{\link{ffsuitable}} }
#! \examples{
#!    message("ffvecapply examples")
#!    x <- ff(vmode="integer", length=1000)
#!    message("loop evaluate expression without returning anything")
#!    ffvecapply(x[i1:i2] <- i1:i2, X=x, VERBOSE=TRUE)
#!    ffvecapply(x[i1:i2] <- i1:i2, X=x, BATCHSIZE=200, VERBOSE=TRUE)
#!    ffvecapply(x[i1:i2] <- i1:i2, X=x, BATCHSIZE=199, VERBOSE=TRUE)
#!    message("lets return the combined expressions as a new ff object")
#!    ffvecapply(i1:i2, N=length(x), VMODE="integer", RETURN=TRUE, BATCHSIZE=200)
#!    message("lets return the combined expressions as a new ram object")
#!    ffvecapply(i1:i2, N=length(x), VMODE="integer", RETURN=TRUE, FF_RETURN=FALSE, BATCHSIZE=200)
#!    message("lets return the combined expressions in existing ff object x")
#!    x[] <- 0L
#!    ffvecapply(i1:i2, N=length(x), VMODE="integer", RETURN=TRUE, FF_RETURN=x, BATCHSIZE=200)
#!    x
#!    message("aggregate and collapse")
#!    ffvecapply(summary(x[i1:i2]), X=x, RETURN=TRUE, CFUN="list", BATCHSIZE=200)
#!    ffvecapply(summary(x[i1:i2]), X=x, RETURN=TRUE, CFUN="crbind", BATCHSIZE=200)
#!    ffvecapply(summary(x[i1:i2]), X=x, RETURN=TRUE, CFUN="cmean", BATCHSIZE=200)
#!
#!    message("how to do colSums with ffrowapply")
#!    x <- ff(1:10000, vmode="integer", dim=c(1000, 10))
#!    ffrowapply(colSums(x[i1:i2,,drop=FALSE]), X=x, RETURN=TRUE, CFUN="list", BATCHSIZE=200)
#!    ffrowapply(colSums(x[i1:i2,,drop=FALSE]), X=x, RETURN=TRUE, CFUN="crbind", BATCHSIZE=200)
#!    ffrowapply(colSums(x[i1:i2,,drop=FALSE]), X=x, RETURN=TRUE, CFUN="csum", BATCHSIZE=200)
#!
#!    message("further ffrowapply examples")
#!    x <- ff(1:10000, vmode="integer", dim=c(1000, 10))
#!    message("loop evaluate expression without returning anything")
#!    ffrowapply(x[i1:i2, ] <- i1:i2, X=x, BATCHSIZE=200)
#!    message("lets return the combined expressions as a new ff object (x unchanged)")
#!    ffrowapply(2*x[i1:i2, ], X=x, RETURN=TRUE, BATCHSIZE=200)
#!    message("lets return a single row aggregate")
#!    ffrowapply(t(apply(x[i1:i2,,drop=FALSE], 1, mean)), X=x, RETURN=TRUE, RETCOL=NULL, BATCHSIZE=200)
#!    message("lets return a 6 column aggregates")
#!    y <- ffrowapply( t(apply(x[i1:i2,,drop=FALSE], 1, summary)), X=x
#!    , RETURN=TRUE, RETCOL=length(summary(0)), BATCHSIZE=200)
#!    colnames(y) <- names(summary(0))
#!    y
#!    message("determine column minima if a complete column does not fit into RAM")
#!    ffrowapply(apply(x[i1:i2,], 2, min), X=x, RETURN=TRUE, CFUN="pmin", BATCHSIZE=200)
#!
#!    message("ffapply examples")
#!    x <- ff(1:720, dim=c(8,9,10))
#!    dimnames(x) <- dummy.dimnames(x)
#!    message("apply function with scalar return value")
#!    apply(X=x[], MARGIN=3:2, FUN=sum)
#!    apply(X=x[], MARGIN=2:3, FUN=sum)
#!    ffapply(X=x, MARGIN=3:2, AFUN="sum", RETURN=TRUE, BATCHSIZE=8)
#!    message("this is what CFUN is based on")
#!    ffapply(X=x, MARGIN=2:3, AFUN="sum", RETURN=TRUE, CFUN="list", BATCHSIZE=8)
#!
#!    message("apply functions with vector or array return value currently have limited support")
#!    apply(X=x[], MARGIN=3:2, FUN=summary)
#!    message("you must use CFUN, the rest is up to you")
#!    y <- ffapply(X=x, MARGIN=3:2, AFUN="summary", RETURN=TRUE, CFUN="list", BATCHSIZE=8)
#!    y
#!    y[[1]]
#!
#!    rm(x); gc()
#! }
#! \keyword{ array }
#! \keyword{ data }


ffvecapply <- function(
  EXPR
, X           = NULL
, N           = NULL
, VMODE       = NULL
, VBYTES      = NULL
, RETURN      = FALSE
, CFUN        = NULL
, USE.NAMES   = TRUE
, FF_RETURN   = TRUE
, BREAK       = '.break'
, FROM        = 'i1'
, TO          = 'i2'
, BATCHSIZE   = .Machine$integer.max
, BATCHBYTES  = getOption("ffbatchbytes")
, VERBOSE     = FALSE
)
{
  if (VERBOSE)
    start.time <- proc.time()
  if (is.null(X)){
    if (is.null(N))
      stop("need N (or X)")
    if (is.null(VBYTES)){
      if (is.null(VMODE))
        stop("need VBYTES (or X)")
      else
        VBYTES <- .rambytes[VMODE]
    }
  }else{
    if (is.null(N))
      N <- length(X)
    if (is.null(VMODE))
      VMODE <- vmode(X)
    if (is.null(VBYTES))
      VBYTES <- .rambytes[VMODE]
  }
  B <- as.integer( min(N, BATCHSIZE, BATCHBYTES %/% VBYTES) )
  bb <- bbatch(N,B)
  if (RETURN){
    nbr <- bb$n + (bb$rb>0)
    i1 <- cumsum(c(1L, rep(bb$b, length.out=nbr-1L)))
    i2 <- cumsum(rep(bb$b, length.out=nbr))
    if (bb$rb) i2[nbr] <- N
    if (is.null(CFUN)){
      if (is.null(VMODE))
        stop("need VMODE (or X) when RETURNing")
      FF_ATTR <- list(vmode = VMODE, length=N, initdata=NULL)
      FF_RET <- ffreturn(FF_RETURN=FF_RETURN, FF_PROTO=X, FF_ATTR=FF_ATTR)
    }else{
      FF_RET <- vector("list", nbr)
      if (USE.NAMES) names(FF_RET) <- paste(i1,i2,sep=":")
    }
  }
  p <- parent.frame()
  assign(FROM, 1L, p)
  assign(TO, bb$b, p)
  e1 <- substitute(i1 <- i2 + 1L, list(i1=as.name(FROM), i2=as.name(TO)))
  e2 <- substitute(i2 <- i2 + b, list(i2=as.name(TO), b=bb$b))
  #slower: two in one, why?? e12 <- substitute({i1 <- i2 + 1L; i2 <- i2 + b}, list(i1=as.name(FROM), i2=as.name(TO), b=b))
  e <- substitute(EXPR)
  for (ib in 1:bb$nb){
    if (VERBOSE){
      cat("ffvecapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
      temp.time <- proc.time()[3]
    }
    temp <- eval(e, p)
    if (RETURN){
      if (is.null(CFUN))
        FF_RET[i1[ib]:i2[ib]] <- temp
      else
        FF_RET[[ib]] <- temp
    }
    if (VERBOSE){
      cat(".. in", proc.time()[3]-temp.time, "seconds\n")
    }
    if (exists(BREAK, p, mode="logical", inherits=FALSE))
      break
    eval(e1, p)
    eval(e2, p)
  }
  if (bb$rb && !exists(BREAK, p, mode="logical", inherits=FALSE)){
    assign(FROM, N-bb$rb+1L, p)
    assign(TO, N, p)
    if (VERBOSE){
      cat("ffvecapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
      temp.time <- proc.time()[3]
    }
    temp <- eval(e, p)
    if (RETURN){
      if (is.null(CFUN))
        FF_RET[i1[nbr]:i2[nbr]] <- temp
      else
        FF_RET[[nbr]] <- temp
    }
    if (VERBOSE){
      cat(".. in", proc.time()[3]-temp.time, "seconds\n")
    }
  }
  rm(list=c(FROM,TO), envir=p)
  if (VERBOSE){
    cat("TOTAL TIME\n")
    print(proc.time() - start.time)
  }
  if (RETURN){
    if (is.null(CFUN) || CFUN=="list")
      FF_RET
    else
      do.call(CFUN, FF_RET)
  }else{
    invisible()
  }
}


ffrowapply <- function(
  EXPR
, X           = NULL
, N           = NULL
, NCOL        = NULL
, VMODE       = NULL
, VBYTES      = NULL
, RETURN      = FALSE
, RETCOL      = NCOL
, CFUN        = NULL
, USE.NAMES   = TRUE
, FF_RETURN   = TRUE
, FROM        = 'i1'
, TO          = 'i2'
, BATCHSIZE   = .Machine$integer.max
, BATCHBYTES  = getOption("ffbatchbytes")
, VERBOSE     = FALSE
)
{
  if (VERBOSE)
    start.time <- proc.time()
  if (is.null(X)){
    if (is.null(N))
      stop("need N (or X)")
    if (is.null(NCOL))
      stop("need NCOL (or X)")
    else
      NCOL <- as.integer(NCOL)
    if (is.null(VBYTES)){
      if (is.null(VMODE))
        stop("need VBYTES (or X)")
      else
        VBYTES <- .rambytes[VMODE]
    }
  }else{
    if (is.null(N))
      N <- nrow(X)
    if (is.null(NCOL))
      NCOL <- ncol(X)
    if (is.null(VMODE))
      VMODE <- vmode(X)
    if (is.null(VBYTES))
      VBYTES <- .rambytes[VMODE]
  }
  B <- as.integer( min(N, BATCHSIZE, BATCHBYTES %/% (VBYTES * NCOL)) )
  if (B<1)
    stop("BATCHBYTES too small")
  bb <- bbatch(N,B)
  if (RETURN){
    nbr <- bb$n + (bb$rb>0)
    i1 <- cumsum(c(1L, rep(bb$b, length.out=nbr-1L)))
    i2 <- cumsum(rep(bb$b, length.out=nbr))
    if (bb$rb) i2[nbr] <- N
    if (is.null(CFUN)){
      if (is.null(VMODE))
        stop("need VMODE (or X) when RETURNing")
      if (is.logical(FF_RETURN) && !FF_RETURN){
        if (is.null(X) || !identical(N, nrow(X)) || !identical(RETCOL, ncol(X))){
          if (is.null(RETCOL)){
            FF_RET <- vector.vmode(VMODE, N)
            names(FF_RET) <- rownames(X)
          }else{
            FF_RET <- vector.vmode(VMODE, N*RETCOL)
            dim(FF_RET)<- c(N,RETCOL)
            rownames(FF_RET) <- rownames(X)
          }
        }else{
          FF_RET <- as.ram(X, vmode=VMODE)
        }
      }else{
        if (is.null(X) || !identical(N, nrow(X)) || !identical(RETCOL, ncol(X))){
          if (is.null(RETCOL)){
            FF_RET <- ff(vmode = VMODE, length=N)
            names(FF_RET) <- rownames(X)
          }else{
            FF_RET <- ff(vmode = VMODE, dim=c(N,RETCOL))
            rownames(FF_RET) <- rownames(X)
          }
        }else{
          FF_ATTR <- list(vmode = VMODE, initdata=NULL)
          FF_RET <- ffreturn(FF_RETURN=FF_RETURN, FF_PROTO=X, FF_ATTR=FF_ATTR)
        }
      }
    }else{
      FF_RET <- vector("list", nbr)
      if (USE.NAMES) names(FF_RET) <- paste(i1,i2,sep=":")
    }
  }
  p <- parent.frame()
  assign(FROM, 1L, p)
  assign(TO, bb$b, p)
  e1 <- substitute(i1 <- i2 + 1L, list(i1=as.name(FROM), i2=as.name(TO)))
  e2 <- substitute(i2 <- i2 + b, list(i2=as.name(TO), b=bb$b))
  #slower: two in one, why?? e12 <- substitute({i1 <- i2 + 1L; i2 <- i2 + b}, list(i1=as.name(FROM), i2=as.name(TO), b=b))
  e <- substitute(EXPR)
  for (ib in 1:bb$nb){
    if (VERBOSE){
      cat("ffrowapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
      temp.time <- proc.time()[3]
    }
    temp <- eval(e, p)
    if (RETURN){
      if (is.null(CFUN)){
        if (is.null(RETCOL))
          FF_RET[i1[ib]:i2[ib]] <- temp
        else
          FF_RET[i1[ib]:i2[ib],] <- temp
      }else
        FF_RET[[ib]] <- temp
    }
    if (VERBOSE){
      cat(".. in", proc.time()[3]-temp.time, "seconds\n")
    }
    eval(e1, p)
    eval(e2, p)
  }
  if (bb$rb){
    assign(FROM, N-bb$rb+1L, p)
    assign(TO, N, p)
    if (VERBOSE){
      cat("ffrowapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
      temp.time <- proc.time()[3]
    }
    temp <- eval(e, p)
    if (RETURN){
      if (is.null(CFUN)){
        if (is.null(RETCOL))
          FF_RET[i1[nbr]:i2[nbr]] <- temp
        else
          FF_RET[i1[nbr]:i2[nbr],] <- temp
      }else
        FF_RET[[nbr]] <- temp
    }
    if (VERBOSE){
      cat(".. in", proc.time()[3]-temp.time, "seconds\n")
    }
  }
  rm(list=c(FROM,TO), envir=p)
  if (VERBOSE){
    cat("TOTAL TIME\n")
    print(proc.time() - start.time)
  }
  if (RETURN){
    if (is.null(CFUN) || CFUN=="list")
      FF_RET
    else
      do.call(CFUN, FF_RET)
  }else{
    invisible()
  }
}


ffcolapply <- function(
  EXPR
, X           = NULL
, N           = NULL
, NROW        = NULL
, VMODE       = NULL
, VBYTES      = NULL
, RETURN      = FALSE
, RETROW      = NROW
, CFUN        = NULL
, USE.NAMES   = TRUE
, FF_RETURN   = TRUE
, FROM        = 'i1'
, TO          = 'i2'
, BATCHSIZE   = .Machine$integer.max
, BATCHBYTES  = getOption("ffbatchbytes")
, VERBOSE     = FALSE
){
  if (VERBOSE)
    start.time <- proc.time()
  if (is.null(X)){
    if (is.null(N))
      stop("need N (or X)")
    if (is.null(NROW))
      stop("need NROW (or X)")
    else
      NROW <- as.integer(NROW)
    if (is.null(VBYTES)){
      if (is.null(VMODE))
        stop("need VBYTES (or X)")
      else
        VBYTES <- .rambytes[VMODE]
    }
  }else{
    if (is.null(N))
      N <- ncol(X)
    if (is.null(NROW))
      NROW <- nrow(X)
    if (is.null(VMODE))
      VMODE <- vmode(X)
    if (is.null(VBYTES))
      VBYTES <- .rambytes[VMODE]
  }
  B <- as.integer( min(N, BATCHSIZE, BATCHBYTES %/% (VBYTES * NROW)) )
  if (B<1)
    stop("BATCHBYTES too small")
  bb <- bbatch(N,B)
  if (RETURN){
    nbr <- bb$n + (bb$rb>0)
    i1 <- cumsum(c(1L, rep(bb$b, length.out=nbr-1L)))
    i2 <- cumsum(rep(bb$b, length.out=nbr))
    if (bb$rb) i2[nbr] <- N
    if (is.null(CFUN)){
      if (is.null(VMODE))
        stop("need VMODE (or X) when RETURNing")
      if (is.logical(FF_RETURN) && !FF_RETURN){
        if (is.null(X) || !identical(N, ncol(X)) || !identical(RETROW, nrow(X))){
          if (is.null(RETROW)){
            FF_RET <- vector.vmode(VMODE, N)
            names(FF_RET) <- colnames(X)
          }else{
            FF_RET <- vector.vmode(VMODE, N*RETROW)
            dim(FF_RET)<- c(RETROW,N)
            colnames(FF_RET) <- colnames(X)
          }
        }else{
          FF_RET <- as.ram(X, vmode=VMODE)
        }
      }else{
        if (is.null(X) || !identical(N, ncol(X)) || !identical(RETROW, nrow(X))){
          if (is.null(RETROW)){
            FF_RET <- ff(vmode = VMODE, length=N)
            names(FF_RET) <- colnames(X)
          }else{
            FF_RET <- ff(vmode = VMODE, dim=c(RETROW, N))
            colnames(FF_RET) <- colnames(X)
          }
        }else{
          FF_ATTR <- list(vmode = VMODE, initdata=NULL)
          FF_RET <- ffreturn(FF_RETURN=FF_RETURN, FF_PROTO=X, FF_ATTR=FF_ATTR)
        }
      }
    }else{
      FF_RET <- vector("list", nbr)
      if (USE.NAMES) names(FF_RET) <- paste(i1,i2,sep=":")
    }
  }
  p <- parent.frame()
  assign(FROM, 1L, p)
  assign(TO, bb$b, p)
  e1 <- substitute(i1 <- i2 + 1L, list(i1=as.name(FROM), i2=as.name(TO)))
  e2 <- substitute(i2 <- i2 + b, list(i2=as.name(TO), b=bb$b))
  #slower: two in one, why?? e12 <- substitute({i1 <- i2 + 1L; i2 <- i2 + b}, list(i1=as.name(FROM), i2=as.name(TO), b=b))
  e <- substitute(EXPR)
  for (ib in 1:bb$nb){
    if (VERBOSE){
      cat("ffcolapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
      temp.time <- proc.time()[3]
    }
    temp <- eval(e, p)
    if (RETURN){
      if (is.null(CFUN)){
        if (is.null(RETROW))
          FF_RET[i1[ib]:i2[ib]] <- temp
        else
          FF_RET[,i1[ib]:i2[ib]] <- temp
      }else
        FF_RET[[ib]] <- temp
    }
    if (VERBOSE){
      cat(".. in", proc.time()[3]-temp.time, "seconds\n")
    }
    eval(e1, p)
    eval(e2, p)
  }
  if (bb$rb){
    assign(FROM, N-bb$rb+1L, p)
    assign(TO, N, p)
    if (VERBOSE){
      cat("ffcolapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
      temp.time <- proc.time()[3]
    }
    temp <- eval(e, p)
    if (RETURN){
      if (is.null(CFUN)){
        if (is.null(RETROW))
          FF_RET[i1[nbr]:i2[nbr]] <- temp
        else
          FF_RET[,i1[nbr]:i2[nbr]] <- temp
      }else
        FF_RET[[nbr]] <- temp
    }
    if (VERBOSE){
      cat(".. in", proc.time()[3]-temp.time, "seconds\n")
    }
  }
  rm(list=c(FROM,TO), envir=p)
  if (VERBOSE){
    cat("TOTAL TIME\n")
    print(proc.time() - start.time)
  }
  if (RETURN){
    if (is.null(CFUN) || CFUN=="list")
      FF_RET
    else
      do.call(CFUN, FF_RET)
  }else{
    invisible()
  }
}


ffapply <- function(
  EXPR        = NULL
, AFUN        = NULL
, MARGIN      = NULL
, X           = NULL
, N           = NULL
, DIM         = NULL
, VMODE       = NULL
, VBYTES      = NULL
, RETURN      = FALSE
, CFUN        = NULL
, USE.NAMES   = TRUE
, FF_RETURN   = TRUE
, IDIM        = 'idim'    # dimension index
, FROM        = 'i1'
, TO          = 'i2'
, BREAK       = '.break'
, BATCHSIZE   = .Machine$integer.max
, BATCHBYTES  = getOption("ffbatchbytes")
, VERBOSE     = FALSE
)
{
  if (VERBOSE)
    start.time <- proc.time()
  if (is.null(X)){
    if (is.null(N)&&is.null(DIM))
      stop("need N or DIM (or X)")
    if (!is.null(N))
      N <- as.integer(N)
    if (!is.null(DIM))
      DIM <- as.integer(DIM)
    if (is.null(VBYTES)){
      if (is.null(VMODE))
        stop("need VBYTES (or X)")
      else
        VBYTES <- .rambytes[VMODE]
    }
  }else{
    if (is.null(N))
      N <- length(X)
    if (is.null(DIM))
      DIM <- dim(X)
    if (is.null(VMODE))
      VMODE <- vmode(X)
    if (is.null(VBYTES))
      VBYTES <- .rambytes[VMODE]
  }


  if (is.null(MARGIN)){
    B <- as.integer( min(N, BATCHSIZE, BATCHBYTES %/% VBYTES) )
    bb <- bbatch(N,B)

    if (RETURN){
      nbr <- bb$nb + (bb$rb>0)
      i1 <- cumsum(c(1L, rep(bb$b, length.out=nbr-1L)))
      i2 <- cumsum(rep(bb$b, length.out=nbr))
      if (bb$rb) i2[nbr] <- N
      if (is.null(CFUN)){
        if (is.null(VMODE))
          stop("need VMODE (or X) when RETURNing")
        FF_ATTR <- list(vmode = VMODE, length=N, initdata=NULL)
        FF_RET <- ffreturn(FF_RETURN=FF_RETURN, FF_PROTO=X, FF_ATTR=FF_ATTR)
      }else{
        FF_RET <- vector("list", nbr)
        if (USE.NAMES) names(FF_RET) <- paste(i1,i2,sep=":")
      }
    }
    p <- parent.frame()
    assign(FROM, 1L, p)
    assign(TO, bb$b, p)
    e1 <- substitute(i1 <- i2 + 1L, list(i1=as.name(FROM), i2=as.name(TO)))
    e2 <- substitute(i2 <- i2 + b, list(i2=as.name(TO), b=bb$b))
   #slower: two in one, why?? e12 <- substitute({i1 <- i2 + 1L; i2 <- i2 + b}, list(i1=as.name(FROM), i2=as.name(TO), b=b))
    if (is.null(AFUN))
      e <- substitute(EXPR)
    else{
      indexargs <- alist(x=)
      names(indexargs) <- NULL
      indexargs[[1]] <- substitute(i1:i2, list(i1=as.symbol(FROM), i2=as.symbol(TO)))
      indexcall <- as.list(call("[", 1))
      indexcall[[2]] <- substitute(X)
      indexcall <- as.call(c(as.list(indexcall), indexargs))
      applycall <- as.list(call("f", 1))
      applycall[[1]] <- substitute(AFUN)
      applycall[[2]]<- indexcall
      e <- as.call(applycall)
    }
    for (ib in 1:bb$nb){
      if (VERBOSE){
        cat("ffapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
        temp.time <- proc.time()[3]
      }
      temp <- eval(e, p)
      if (RETURN){
        if (is.null(CFUN))
          FF_RET[i1[ib]:i2[ib]] <- temp
        else
          FF_RET[[ib]] <- temp
      }
      if (VERBOSE){
        cat(".. in", proc.time()[3]-temp.time, "seconds\n")
      }
      if (exists(BREAK, p, mode="logical", inherits=FALSE))
        break
      eval(e1, p)
      eval(e2, p)
    }
    if (bb$rb && !exists(BREAK, p, mode="logical", inherits=FALSE)){
      assign(FROM, N-bb$rb+1L, p)
      assign(TO, N, p)
      if (VERBOSE){
        cat("ffapply", paste(get(FROM,p),":",get(TO,p),"{", N, "}", sep="", collapse=", "), "..")
        temp.time <- proc.time()[3]
      }
      temp <- eval(e, p)
      if (RETURN){
        if (is.null(CFUN))
          FF_RET[i1[nbr]:i2[nbr]] <- temp
        else
          FF_RET[[nbr]] <- temp
      }
      if (VERBOSE){
        cat(".. in", proc.time()[3]-temp.time, "seconds\n")
      }
    }
    rm(list=c(FROM,TO), envir=p)
    if (VERBOSE){
      cat("TOTAL TIME\n")
      print(proc.time() - start.time)
    }
    if (RETURN){
      if (is.null(CFUN) || CFUN=="list")
        FF_RET
      else
        do.call(CFUN, FF_RET)
    }else{
      invisible()
    }
  }else{

    if (is.null(DIM))
      stop("specified MARGIN but is.null(DIM)")
    if (length(DIM)<length(MARGIN))
      stop("DIM is shorter than MARGIN")
    MARGIN <- as.integer(MARGIN)
    ndim <- length(DIM)
    mdim <- DIM[MARGIN]
    nmarg <- length(MARGIN)
    MARGIN.rem <- (1:ndim)[-MARGIN]
    dims <- c(rev(MARGIN), rev(MARGIN.rem)) # slowest rotating first
    chunksize <- rev(cumprod(c(1, rev(DIM[dims])))) # cumprod from last to first
    N <- chunksize[1]
    chunksize <- chunksize[-1] # remove n at the beginning, now we have chunksize at each dimension
    bytes <- .rambytes[VMODE]
    remsize <- chunksize[nmarg]
    if (remsize>BATCHSIZE || (remsize*bytes)>BATCHBYTES)
      warning("ignoring that finest batch size|bytes ", remsize, "|", remsize*bytes, " > BATCHSIZE|BATCHBYTES ", BATCHSIZE, "|", BATCHBYTES)
    margsize <- chunksize[nmarg:1L]

    B <- as.integer(pmax(1,pmin(mdim, BATCHSIZE %/% margsize, BATCHBYTES %/% (margsize * bytes)))) # batch sizes per dimenion
    bb <- bbatch(mdim,B)
    b <- bb$b
    nb <- bb$nb
    rb <- bb$rb
    if (RETURN){
      nbr <- nb + (rb>0)
      if (is.null(CFUN)){
        if (is.null(VMODE))
          stop("need VMODE (or X) when RETURNing")
        FF_ATTR <- list(vmode=VMODE, dim=mdim, dimnames=dimnames(X)[MARGIN], dimorder=order(match(mdim, dimorder(X))), initdata=NULL)
        FF_RET <- ffreturn(FF_RETURN=FF_RETURN, FF_PROTO=X, FF_ATTR=FF_ATTR)
      }else{
        FF_RET <- array(list(), nbr)
        if (USE.NAMES){
          i12 <- lapply(1:nmarg, function(i){
              i1 <- cumsum(c(1L, rep(b[i], length.out=nbr[i]-1L)))
              i2 <- cumsum(rep(b[i], length.out=nbr[i]))
              if (rb[i]) i2[nbr[i]] <- mdim[i]
              list(i1=i1,i2=i2)
            })
          dimnames(FF_RET) <- lapply(i12, function(i)paste(i$i1,i$i2,sep=":"))
        }
      }
    }

    if (is.null(AFUN))
      applycall <- substitute(EXPR)
    else{
      aarg <- alist(x=)
      names(aarg) <- NULL
      indexargs <- rep(aarg,ndim)
      for (i in 1:nmarg)
        indexargs[[MARGIN[i]]] <- substitute(i1[i]:i2[i], list(i=i, i1=as.symbol(FROM), i2=as.symbol(TO)))
      indexargs <- c(indexargs, alist(drop=FALSE))
      indexcall <- as.list(call("[", 1))
      indexcall[[2]] <- substitute(X)
      indexcall <- as.call(c(as.list(indexcall), indexargs))
      applycall <- as.list(call("apply", X=1, MARGIN=MARGIN, FUN=AFUN))
      applycall$X <- indexcall
      applycall <- as.call(applycall)
    }
    if (VERBOSE){
      e <- substitute({
        cat("ffapply", paste(i1,":",i2,"{", mdim, "}", sep="", collapse=", "), "..")
        1L
      }, list(IDIM=as.name(IDIM), i1=as.name(FROM), i2=as.name(TO), mdim=mdim))
      e[[3]] <- applycall
    }else{
      e <- applycall
    }
    e12 <- substitute({
      .ffapply.dimexhausted <- i2[IDIM]>=mdim[IDIM]
      while(.ffapply.dimexhausted){
        if (IDIM<nmarg){
          i1[IDIM] <- 1L
          i2[IDIM] <- b[IDIM]
          IDIM <- IDIM + 1
          .ffapply.dimexhausted <- i2[IDIM]>=mdim[IDIM]
        }else{
          .ffapply.dimexhausted <- FALSE
          BREAK <- TRUE
        }
      }
      if (!exists(BREAK, mode="logical", inherits=FALSE)){
        i1[IDIM] <- i2[IDIM] + 1L
        i2[IDIM] <- min(mdim[IDIM], i2[IDIM] + b[IDIM])
        IDIM <- 1L
        1L
      }else{
        NULL
      }
    }, list(IDIM=as.name(IDIM), i1=as.name(FROM), i2=as.name(TO), nmarg=nmarg, MARGIN=MARGIN, mdim=mdim, b=b, VERBOSE=VERBOSE, BREAK=BREAK))
    e12[[4]][[3]][[5]] <- e
    ib <- 1L
    p <- parent.frame()
    assign(IDIM, 1L, p)
    i1 <- rep(1L, nmarg); i1[1] <- 0L
    i2 <- b; i2[1] <- 0L
    assign(FROM, i1, p)
    assign(TO, i2, p)
    goon <- TRUE
    if (VERBOSE){
      temp.time <- proc.time()[3]
      while(goon){
        temp <- eval(e12, p)
        goon <- !exists(BREAK, p, mode="logical", inherits=FALSE)
        if (goon){
          temp2.time <- proc.time()[3]
          cat(".. in", temp2.time-temp.time, "seconds\n")
          temp.time <- temp2.time
          if (RETURN){
            if (is.null(CFUN)){
              i1 <- get(FROM, p)
              i2 <- get(TO, p)
              FF_RET <- do.call( "[<-", c( list(FF_RET), lapply(1:nmarg, function(i)substitute(e1:e2, list(e1=i1[i], e2=i2[i]))), list(value=temp) ) )
            }else
              FF_RET[[ib]] <- temp
            ib <- ib + 1L
          }
        }
      }
    }else{
      while(goon){
        temp <- eval(e12, p)
        goon <- !exists(BREAK, p, mode="logical", inherits=FALSE)
        if (goon){
          if (RETURN){
            if (is.null(CFUN)){
              i1 <- get(FROM, p)
              i2 <- get(TO, p)
              FF_RET <- do.call( "[<-", c( list(FF_RET), lapply(1:nmarg, function(i)substitute(e1:e2, list(e1=i1[i], e2=i2[i]))), list(value=temp) ) )
            }else
              FF_RET[[ib]] <- temp
            ib <- ib + 1L
          }
        }
      }
    }
    rm(list=c(FROM,TO,".ffapply.dimexhausted", IDIM), envir=p)
    if (exists(BREAK, p, mode="logical", inherits=FALSE))
      rm(list=BREAK, envir=p)
    if (VERBOSE){
      cat("TOTAL TIME\n")
      print(proc.time() - start.time)
    }
    if (RETURN){
      if (is.null(CFUN) || CFUN=="list")
        FF_RET
      else
        do.call(CFUN, FF_RET)
    }else{
      invisible()
    }

  }

}


if (FALSE){
  library(ff)

  d <- c(4,4,4,4)
  a <- array(1:(cumprod(d)[length(d)]), d)
  dimnames(a) <- lapply(1:length(dim(a)), function(i)paste(letters[i], 1:(dim(a)[i]), sep=""))
  apply(a, 1:2, paste, collapse="|")
  apply(a, 1:2, sum)
  ffapply(a, apply(a[i1[1]:i2[1],i1[2]:i2[2],,,drop=FALSE], 3:4, sum), 1:2, BATCHSIZE=16, VERBOSE=T, RETURN="list")
  ffapply(a, apply(a[i1[1]:i2[1],i1[2]:i2[2],,,drop=FALSE], 3:4, sum), 1:2, BATCHSIZE=16, VERBOSE=T, RETURN=csum)
  ffapply(a, apply(a[i1[1]:i2[1],i1[2]:i2[2],,,drop=FALSE], 3:4, sum), 1:2, BATCHSIZE=16, VERBOSE=T, RETURN=c)

  ffapply(a, apply(a[,,i1[1]:i2[1],i1[2]:i2[2],drop=FALSE], 1:2, sum), 3:4, BATCHSIZE=16, VERBOSE=T, RETURN=csum)

  ffapply(a, apply(a[i1[1]:i2[1],i1[2]:i2[2],,,drop=FALSE], 3:4, sum), 1:2, BATCHSIZE=16, VERBOSE=T, RETURN="list")

  n <- 1000000
  x <- ff(1, dim=c(n,2))
  y <- ffapply(x, rowSums(x[i1:i2,,drop=FALSE]), MARGIN=1, RETURN="unlist", BATCHSIZE=10000)
  z <- ffrowapply(x, rowSums(x[i1:i2,,drop=FALSE]), RETURN="unlist", BATCHSIZE=10000)

  system.time(ffapply(x, i2, MARGIN=1, RETURN="unlist", BATCHSIZE=1000))
  system.time(ffrowapply(x, i2, RETURN="unlist", BATCHSIZE=1000))

  n <- 1000000
  x <- ff(0:0, n)
  system.time(ffvecapply(x[i1:i2], x, BATCHSIZE=1000))

  x <- 1:100
  print(ffapply(x,{x[i1:i2] <- x[i1:i2] + 1;c(i1,i2)}, BATCHSIZE=10, RETURN=csum))
  x[1:10]
}
