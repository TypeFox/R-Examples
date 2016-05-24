# Chunking utilities for bit and ff
# (c) 2007-2009 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-09-03
# Last changed: 2007-10-25

# source("D:/mwp/eanalysis/bit/R/chunkutil.R")



#! \name{bbatch}
#! \alias{bbatch}
#! \title{ Balanced Batch sizes }
#! \description{
#!   \command{bbatch} calculates batch sizes so that they have rather balanced sizes than very different sizes
#! }
#! \usage{
#! bbatch(N, B)
#! }
#! \arguments{
#!   \item{N}{ total size }
#!   \item{B}{ desired batch size }
#! }
#! \value{
#!   a list with components
#!   \item{ b }{ the batch size }
#!   \item{ nb }{ the number of batches }
#!   \item{ rb }{ the size of the rest }
#! }
#! \details{
#!   Tries to have \code{rb==0} or \code{rb} as close to \code{b} as possible while guaranteing that \code{rb < b && (b - rb) <= min(nb, b)}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{repfromto}}, \code{\link[ff]{ffvecapply}} }
#! \examples{
#!   bbatch(100, 24)
#! }
#! \keyword{ IO }
#! \keyword{ data }

# non-vectorized
#bbatch <-
#function(N,B){
#  N <- as.integer(N)
#  B <- as.integer(B)
#  RB <- N %% B
#  NB <- N %/% B
#  if (RB){
#    cc <- min((B - RB) %/% NB, (B - RB) %/% (NB + 1L))
#    if (cc){
#      rb <- RB + cc * NB
#      b <- B - cc
#      if (rb==b){
#        return(list(b=b, nb=NB+1L, rb=0L))
#      }else{
#        return(list(b=b, nb=NB, rb=rb))
#      }
#    }else{
#      return(list(b=B, nb=NB, rb=RB))
#    }
#  }else{
#    return(list(b=B, nb=NB, rb=RB))
#  }
#}

# balanced batchsizes
bbatch <- function(N,B){
  if (any(B<1))
    stop("B too small")
  N <- as.integer(N)
  B <- as.integer(B)
  RB <- N %% B
  NB <- N %/% B
  cc <- pmin((B - RB) %/% NB, (B - RB) %/% (NB + 1L))
  cc[RB==0 | NB == 0] <- 0
  i <- cc > 0
  RB[i] <- RB[i] + cc[i] * NB[i]
  B[i] <- B[i] - cc[i]
  j <- i & (RB[i] == B[i])
  NB[i & j] <- NB[i & j] + 1L
  RB[i & j] <- 0L
  i <- (RB>0) & (NB == 0)
  B[i] <- RB[i]
  NB[i] <- 1L
  RB[i] <- 0L
  return(list(b=B, nb=NB, rb=RB))
}

#! \name{repfromto}
#! \alias{repfromto}
#! \alias{repfromto<-}
#! \title{ Virtual recycling }
#! \description{
#!   \command{repfromto} virtually recylcles object \code{x} and cuts out positions \code{from .. to}
#! }
#! \usage{
#! repfromto(x, from, to)
#! repfromto(x, from, to) <- value
#! }
#! \arguments{
#!   \item{x}{ an object from which to recycle }
#!   \item{from}{ first position to return }
#!   \item{to}{ last position to return }
#!   \item{value}{ value to assign }
#! }
#! \details{
#!   \code{repfromto} is a generalization of \code{\link{rep}}, where \code{rep(x, n) == repfromto(x, 1, n)}.
#!   You can see this as an R-side (vector) solution of the \code{mod_iterate} macro in arithmetic.c
#! }
#! \value{
#!   a vector of length \code{from - to + 1}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{rep}}, \code{\link[ff]{ffvecapply}} }
#! \examples{
#!   message("a simple example")
#!   repfromto(0:9, 11, 20)
#! }
#! \keyword{ IO }
#! \keyword{ data }

repfromto <- function(x, from, to){
  nx <- length(x)
  if (nx){
    from <- as.integer(from)
    to <- as.integer(to)
    if (to>nx){
      N <- to - from + 1L
      from <- (from-1L)%%nx + 1L
      to <- to%%nx
      # NOTE: fetch in sequence pre-main-post in case is.ff(x)
      if (from<=to && N<nx){
        ret <- x[from:to]
      }else{
        pre <- x[from:nx]
        nrep <- (N - length(pre) - to) %/%nx
        main <- if (nrep) rep(x[1:nx], nrep) else NULL
        post <- if (to) x[1:to] else NULL
        ret <- c(pre, main, post)
      }
    }else{
      ret <- x[from:to]
    }
    a <- attributes(x[1])
    a$names <- NULL
    attributes(ret) <- a
    return(ret)
  }else{
    return(NA[from:to])
  }
}

"repfromto<-" <- function(x, from, to, value){
  x[from:to] <- value
  x
}


if (FALSE){
  x <- 1:10
  for (n in 1:20)
  for (i1 in 1:30){
    i2 <- i1+n-1
    cat(i1,i2,"|",repfromto(x,i1,i2), "\n")
  }
}

if (FALSE){
  intseq <- function(from=NULL, to=NULL, by=NULL, length.out=NULL, along.with=NULL){
    if (is.null(from)){
      if (is.null(to))
        stop("need 'from' or 'to'")
      else
        to <- as.integer(to)
      if (is.null(by))
        by <- 1L
      else
        by <- as.integer(by)
    }else{
      from <- as.integer(from)
      if (is.null(to)){
        if (is.null(by))
          by <- 1L
        else
          by <- as.integer(by)
      }else{
        to <- as.integer(to)
        n <- to - from
        if (is.null(by)){
          if (to<from)
            by = -1L
          else
            by = 1L
        }else{
          by <- as.integer(by)
          if (n){
            if (sign(n) != sign(by))
              stop("wrong sign of by")
          }else
            return(from)  # to == from
        }
      }
    }

    if (is.null(length.out)){
      if (is.null(along.with)){
        if (is.null(to) || is.null(from))
          stop("not enough info to guess the length.out")
        else{
          length.out <- n %/% by + 1L
        }
      }else{
        length.out <- length(along.with)
      }
    }else{
      length.out <- as.integer(length.out)
    }
    if (length.out){
      if (length.out==1L)
        from
      else
        cumsum(c(from, rep(by, length.out-1L)))
    }else
      integer()
  }
}


#! \name{chunk}
#! \alias{chunk}
#! \alias{chunk.default}
#! \title{ Chunked range index }
#! \description{
#!   creates a sequence of range indexes using a syntax not completely unlike 'seq'
#! }
#! \usage{
#! chunk(\dots)
#! \method{chunk}{default}(from = NULL, to = NULL, by = NULL, length.out = NULL, along.with = NULL
#! , overlap = 0L, method = c("bbatch", "seq"), maxindex = NA, \dots)
#! }
#! \arguments{
#!   \item{from}{ the starting value of the sequence. }
#!   \item{to}{ the (maximal) end value of the sequence. }
#!   \item{by}{ increment of the sequence }
#!   \item{length.out}{ desired length of the sequence. }
#!   \item{along.with}{ take the length from the length of this argument. }
#!   \item{overlap}{ number of values to overlap (will lower the starting value of the sequence, the first range becomes smaller }
#!   \item{method}{ default 'bbatch' will try to balance the chunk size, see \code{\link{bbatch}}, 'seq' will create chunks like \code{\link[base]{seq}} }
#!   \item{maxindex}{ passed to \code{\link{ri}} }
#!   \item{\dots}{ ignored }
#! }
#! \details{
#!   \code{chunk} is generic, the default method is described here, other methods that automatically consider RAM needs are provided with package 'ff', see for example \code{\link[ff]{chunk.ffdf}}
#! }
#! \section{available methods}{
#!   \code{chunk.default}, \code{\link[ff]{chunk.bit}}, \code{\link[ff]{chunk.ff_vector}}, \code{\link[ff]{chunk.ffdf}}
#! }
#! \value{
#!   \code{chunk.default} returns a list of \code{\link{ri}} objects representing chunks of subscripts
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ri}},  \code{\link[base]{seq}}, \code{\link{bbatch}} }
#! \examples{
#!   chunk(1, 100, by=30)
#!   chunk(1, 100, by=30, method="seq")
#!    \dontrun{
#! require(foreach)
#! m <- 10000
#! k <- 1000
#! n <- m*k
#! message("Four ways to loop from 1 to n. Slowest foreach to fastest chunk is 1700:1 
#! on a dual core notebook with 3GB RAM\n")
#! z <- 0L; 
#! print(k*system.time({it <- icount(m); foreach (i = it) \%do\% { z <- i; NULL }}))
#! z
#!
#! z <- 0L
#! print(system.time({i <- 0L; while (i<n) {i <- i + 1L; z <- i}}))
#! z
#!
#! z <- 0L
#! print(system.time(for (i in 1:n) z <- i))
#! z
#!
#! z <- 0L; n <- m*k; 
#! print(system.time(for (ch in chunk(1, n, by=m)){for (i in ch[1]:ch[2])z <- i}))
#! z
#!
#! message("Seven ways to calculate sum(1:n). 
#!  Slowest foreach to fastest chunk is 61000:1 on a dual core notebook with 3GB RAM\n")
#! print(k*system.time({it <- icount(m); foreach (i = it, .combine="+") \%do\% { i }}))
#! 
#! z <- 0; 
#! print(k*system.time({it <- icount(m); foreach (i = it) \%do\% { z <- z + i; NULL }}))
#! z
#!
#! z <- 0; print(system.time({i <- 0L;while (i<n) {i <- i + 1L; z <- z + i}})); z
#!
#! z <- 0; print(system.time(for (i in 1:n) z <- z + i)); z
#!
#! print(system.time(sum(as.double(1:n))))
#!
#! z <- 0; n <- m*k
#! print(system.time(for (ch in chunk(1, n, by=m)){for (i in ch[1]:ch[2])z <- z + i}))
#! z
#!
#! z <- 0; n <- m*k
#! print(system.time(for (ch in chunk(1, n, by=m)){z <- z+sum(as.double(ch[1]:ch[2]))}))
#! z
#!    }
#! }
#! \keyword{ data }

chunk <- function(...)
  UseMethod("chunk")

chunk.default <- function(
  from = NULL
, to = NULL
, by = NULL
, length.out = NULL
, along.with = NULL
, overlap = 0L
, method=c("bbatch","seq")
, maxindex = NA
, ...
)
{
  method <- match.arg(method)
  if (!is.null(along.with)){
    if (is.null(from))
      from <- 1L
    else{

      if (length(from)==1)
        from <- as.integer(from)
      else
        stop("'from' must be scalar")
    }
    if (is.null(to))
      to <- length(along.with)
    else{
      if (length(to)==1)
        to <- as.integer(to)
      else
        stop("'to' must be scalar")
    }
  }
  if (length(from)==1)
    from <- as.integer(from)
  else
    stop("'from' must be scalar")
  if (length(to)==1)
    to <- as.integer(to)
  else
    stop("'to' must be scalar")

  if (to<from)
    stop("to < from")
  N <- to - from + 1L

  if (is.null(by)){
    if (is.null(length.out))
      stop("need either 'by' or 'length.out'")
    else{
      if (length(length.out)==1){
        length.out <- as.integer(length.out)
        if (length.out>N)
          length.out <- N
        by <- N %/% length.out
      }else
        stop("'length.out' must be scalar")
    }
  }else{
    if (length(by)==1){
      by <- as.integer(by)
      if (by<1)
        stop("'by' must be > 0")
      length.out <- (N - 1L) %/% by + 1L
    }else
      stop("'by' must be scalar")
  }

  if (method=="bbatch")
    by <- bbatch(N, by)$b

  if (length.out>1L){
    from <- cumsum(c(from, rep(by, length.out - 1L)))
    to <- c(from[-1], from[1] + N) - 1L  # fixed by Edwin de Jonge, 18.1.2011
    if (overlap>0)
      from[-1] <- from[-1] - overlap
  }
  n <- length(from)
  s <- seq_len(n)
  ret <- vector("list", n)
  for (i in s){
    ret[[i]] <- ri(from[i], to[i], maxindex)
  }
  ret
}



#! \name{vecseq}
#! \alias{vecseq}
#! \title{ Vectorized Sequences }
#! \description{
#!   \command{vecseq} returns concatenated multiple sequences
#! }
#! \usage{
#!  vecseq(x, y=NULL, concat=TRUE, eval=TRUE)
#! }
#! \arguments{
#!   \item{x}{ vector of sequence start points }
#!   \item{y}{ vector of sequence end points (if \code{is.null(y)} then \code{x} are taken as endpoints, all starting at 1) }
#!   \item{concat}{ vector of sequence end points (if \code{is.null(y)} then \code{x} are taken as endpoints, all starting at 1) }
#!   \item{eval}{ vector of sequence end points (if \code{is.null(y)} then \code{x} are taken as endpoints, all starting at 1) }
#! }
#! \details{
#!   This is a generalization of \code{\link{sequence}} in that you can choose sequence starts other than 1 and also have options to no concat and/or return a call instead of the evaluated sequence.
#! }
#! \value{
#!   if \code{concat==FALSE} and \code{eval==FALSE} a list with n calls that generate sequences \cr
#!   if \code{concat==FALSE} and \code{eval==TRUE } a list with n sequences \cr
#!   if \code{concat==TRUE } and \code{eval==FALSE} a single call generating the concatenated sequences \cr
#!   if \code{concat==TRUE } and \code{eval==TRUE } an integer vector of concatentated sequences
#! }
#! \author{ Angelo Canty, Jens Oehlschlägel }
#! \seealso{ \code{\link{:}}, \code{\link{seq}}, \code{\link{sequence}} }
#! \examples{
#!   sequence(c(3,4))
#!   vecseq(c(3,4))
#!   vecseq(c(1,11), c(5, 15))
#!   vecseq(c(1,11), c(5, 15), concat=FALSE, eval=FALSE)
#!   vecseq(c(1,11), c(5, 15), concat=FALSE, eval=TRUE)
#!   vecseq(c(1,11), c(5, 15), concat=TRUE, eval=FALSE)
#!   vecseq(c(1,11), c(5, 15), concat=TRUE, eval=TRUE)
#! }
#! \keyword{ manip }

vecseq <- function(x, y=NULL, concat=TRUE, eval=TRUE){
        if (missing(y)){
          y <- x
          x <- 1L
				}
        if (concat){
          if (eval){
						# pure R version was: eval(parse(text=paste("c(",paste(x,y,sep=":",collapse=","),")")))
						# now calling C-code
						nx <- length(x)
						ny <- length(y)
						if (nx<ny)
							x <- rep(as.integer(x), length.out=ny)
						if (ny<nx)
							y <- rep(as.integer(y), length.out=nx)
            .Call("R_bit_vecseq", as.integer(x),  as.integer(y), PACKAGE="bit") 
          }else
            parse(text=paste("c(",paste(x,y,sep=":",collapse=","),")"))[[1]]
        }else{
          if (eval)
            eval(parse(text=paste("list(",paste(x,y,sep=":",collapse=","),")")))
          else
            as.list(parse(text=paste(x,y,sep=":")))
        }
}
