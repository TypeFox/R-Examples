#-- refdata.r ------------------------------------
# Jens Oehlschlaegel
# created:            28.09.03
# performance tested: 28.09.03
# documented:         28.09.03
# changes 06.03.2004: optimal.index(oi) now is idempotent, also for oi=NULL which is interpreted as 'no indexing required'
# changes 06.03.2004: new function need.index(oi) returning TRUE if indexing is required (returning FALSE for NULL)
# changes 06.03.2004: new function posi.index(oi) converting optimal index into positive integers (does not make sense for NULL and non-optimal indices)
# (gpl) 2003
#---------------------------------------------

#source("d:/MWP/eAnalysis/ref/R/refdata.r")

#! \name{optimal.index}
#! \alias{optimal.index}
#! \alias{need.index}
#! \alias{posi.index}
#! \title{ creating standardized, memory optimized index for subsetting }
#! \description{
#!   Function \code{optimal.index} converts an index specification of type \{logical, integer, -integer, character\} into one of \{integer, -integer\} whatever is smaller.
#!   Function \code{need.index} returns TRUE if the index does represent a subset (and thus indexing is needed).
#!   Function \code{posi.index} returns positive integers representing the (sub)set.
#! }
#! \usage{
#! optimal.index(i, n=length(i.names), i.names = names(i), i.previous = NULL, strict = TRUE)
#! need.index(oi)
#! posi.index(oi)
#! }
#! \arguments{
#!   \item{i}{ the original one-dimensional index }
#!   \item{n}{ length of the indexed dimension  (potential iMax if i where integer), not necessary if names= is given }
#!   \item{i.names}{ if i is character then names= represents the names of the indexed dimension }
#!   \item{i.previous}{ if i.previous= is given, the returned index represents \code{x[i.previous][i] == x[optimal.index]} rather than \code{x[i] == x[optimal.index]} }
#!   \item{strict}{ set to FALSE to allow for NAs and duplicated index values, but see details }
#!   \item{oi}{ a return value of \code{optimal.index} }
#! }
#! \details{
#!   When strict=TRUE it is expected that i does not contain NAs and no duplicated index values. Then \code{ identical(x[i], x[optimal.index(i, n=length(x), i.names=names(x))$i]) == TRUE } . \cr
#!   When strict=FALSE i may contain NAs and/or duplicated index values. In this case length optimisation is not performed and optimal.index always returns positive integers.
#! }
#! \note{
#!   \code{need.index(NULL)} is defined and returns FALSE. This allows a function to have an optional parameter oi=NULL and to determine the need of subsetting in one reqest.
#! }
#! \value{
#!   \code{optimal.index} returns the index oi with attributes n=n and ni=length(x[optimal.index]) (which is n-length(i) when i is negative).
#!   \code{need.index} returns a logical scalar
#!   \code{posi.index}  returns a vector of positive integers (or integer(0))
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{refdata}}
#!           \cr please ignore the following unpublished links: ids2index, shift.index, startstop2index
#! }
#! \examples{
#!   l <- letters
#!   names(l) <- letters
#!   stopifnot({i <- 1:3 ; identical(l[i], l[optimal.index(i, n=length(l))])})
#!   stopifnot({i <- -(4:26) ; identical(l[i], l[optimal.index(i, n=length(l))])})
#!   stopifnot({i <- c(rep(TRUE, 3), rep(FALSE, 23)) 
#!     identical(l[i], l[optimal.index(i, n=length(l))])})
#!   stopifnot({i <- c("a", "b", "c")
#!     identical(l[i], l[optimal.index(i, i.names=names(l))])})
#!   old.options <- options(show.error.messages=FALSE)
#!     stopifnot(inherits(try(optimal.index(c(1:3, 3), n=length(l))), "try-error"))
#!   options(old.options)
#!   stopifnot({i <- c(1:3, 3, NA)
#!     identical(l[i], l[optimal.index(i, n=length(l), strict=FALSE)])})
#!   stopifnot({i <- c(-(4:26), -26)
#!     identical(l[i], l[optimal.index(i, n=length(l), strict=FALSE)])})
#!   stopifnot({i <- c(rep(TRUE, 3), rep(FALSE, 23), TRUE, FALSE, NA)
#!     identical(l[i], l[optimal.index(i, n=length(l), strict=FALSE)])})
#!   stopifnot({i <- c("a", "b", "c", "a", NA)
#!     identical(l[i], l[optimal.index(i, i.names=names(l), strict=FALSE)])})
#!   rm(l)
#! }
#! \keyword{ utilities }
#! \keyword{ manip }


optimal.index <- function(i, n=length(i.names), i.names=names(i), i.previous=NULL, strict=TRUE){
  # idempotence
  if (is.null(i) || (!is.null(attr(i, "ni")) && !is.null(attr(i, "n"))))
    return(i)
  index <- seq(length=n)
  if (!is.null(i.names) && (is.character(i) || is.character(i.previous))){
    names(index) <- i.names
  }
  if (is.null(i.previous)){
    i <- as.vector(index[i])
  }else{
    i <- as.vector(index[i.previous][i])
  }
  ni <- length(i)
  if (strict){
    if (any(is.na(i)))
      stop("not all index values matched")
    if (any(duplicated(i)))
      stop("index not unique, you can only select subsets, not supersets")
    if ( ni > (n/2) ){
      # return either positive or negative indices, whatever is smaller
      i <- -((seq(length=n))[is.na(match(seq(length=n), i))])
    }
  }
  attributes(i) <- list(n=n, ni=ni)
  i
}

posi.index <- function(oi){
  if (attr(oi, "ni")){
    if (length(oi)){
      if (oi[1]>0)
        oi
      else
        (1:attr(oi, "n"))[oi]
    }else{
      1:attr(oi, "n")
    }
  }else{
    integer(0)
  }
}

need.index <- function(oi)
{
  if (length(oi))
    return(TRUE)
  ni <- attr(oi, "ni")
  !is.null(ni) && ni==0
}




#! \name{refdata}
#! \alias{refdata}
#! \alias{derefdata}
#! \alias{derefdata<-}
#! \alias{[.refdata}
#! \alias{[<-.refdata}
#! \alias{[[.refdata}
#! \alias{[[<-.refdata}
#! \alias{$.refdata}
#! \alias{$<-.refdata}
#! \alias{dim.refdata}
#! \alias{dimnames.refdata}
#! \alias{row.names.refdata}
#! \alias{names.refdata}
#! \alias{print.refdata}
#! \title{ subsettable reference to matrix or data.frame }
#! \description{
#!   Function \code{refdata} creates objects of class refdata which behave not totally unlike matrices or data.frames but allow for much more memory efficient handling.
#! }
#! \usage{
#! # -- usage for R CMD CHECK, see below for human readable version -----------
#! refdata(x)
#! derefdata(x)
#! derefdata(x) <- value
#!  \method{[}{refdata}(x, i = NULL, j = NULL, drop = FALSE, ref = FALSE)
#!  \method{[}{refdata}(x, i = NULL, j = NULL, ref = FALSE) <- value
#!  \method{dim}{refdata}(x)
#!  \method{dimnames}{refdata}(x)
#!  \method{row.names}{refdata}(x)
#!  \method{names}{refdata}(x)
#!
#! # -- most important usage for human beings  --------------------------------
#! # rd <- refdata(x)             # create reference
#! # derefdata(rd)                # retrieve original data
#! # derefdata(rd) <- value       # modify original data
#! # rd[]                         # get all (current) data
#! # rd[i, j]                     # get part of data
#! # rd[i, j, ref=TRUE]           # get new reference on part of data
#! # rd[i, j]           <- value  # modify / create local copy
#! # rd[i, j, ref=TRUE] <- value  # modify original data (respecting subsetting history)
#! # dim(rd)                      # dim of (subsetted) data
#! # dimnames(rd)                 # dimnames of (subsetted) data
#! }
#! \arguments{
#!   \item{x}{ a matrix or data.frame or any other 2-dimensional object that has operators "[" and "[<-" defined }
#!   \item{i}{ row index }
#!   \item{j}{ col index }
#!   \item{ref}{ FALSE by default. In subsetting: FALSE returns data, TRUE returns new refdata object. In assignments: FALSE modifies a local copy and returns a refdata object embedding it, TRUE modifies the original. }
#!   \item{drop}{ FALSE by default, i.e. returned data have always a dimension attribute. TRUE drops dimension in some cases, the exact result depends on whether a \code{\link{matrix}} or \code{\link{data.frame}} is embedded }
#!   \item{value}{ some value to be assigned }
#! }
#! \details{
#!   Refdata objects store 2D-data in one environment and index information in another environment. Derived refdata objects usually share the data environment but not the index environment. \cr
#!   The index information is stored in a standardized and memory efficient form generated by \code{\link{optimal.index}}. \cr
#!   Thus refdata objects can be copied and subsetted and even modified without duplicating the data in memory. \cr
#!   Empty square bracket subsetting (\code{rd[]}) returns the data, square bracket subsetting (\code{rd[i, j]}) returns subsets of the data as expected. \cr
#!   An additional argument (\code{rd[i, j, ref=TRUE]}) allows to get a reference that stores the subsetting indices. Such a reference behaves transparently as if a smaller matrix/data.frame would be stored and can be subsetted again recursively.
#!   With ref=TRUE indices are always interpreted as row/col indices, i.e. \code{x[i]} and \code{x[cbind(i, j)]} are undefined (and raise stop errors) \cr
#!   Standard square bracket assignment (\code{rd[i, j] <- value}) creates a reference to a locally modified copy of the (potentially subsetted) data. \cr
#!   An additional argument (\code{rd[i, j, ref=TRUE] <- value}) allows to modify the original data, properly recognizing the subsetting history. \cr
#!   A method \code{\link{dim}(refdata)} returns the dim of the (indexed) data. \cr
#!   A \code{\link{dimnames}(refdata)} returns the dimnames of the (indexed) data. \cr
#! }
#! \note{
#!   The refdata code is currently R only (not implemented for S+). \cr
#!   Please note the following differences to matrices and dataframes: \cr
#!   \describe{
#!      \item{\code{x[]}}{you need to write \code{x[]} instead of \code{x} in order to get all current data}
#!      \item{\code{drop=FALSE}}{by default drop=FALSE which gives consistent behaviour for matrices and data.frames. You can use the $- or [[-operator to extract single column vectors which are granted to be of a consistent data type. However, currently $ and [[ are only wrappers to [. They might be performance tuned in later versions.}
#!      \item{\code{x[i]}}{single index subsetting is not defined, use \code{x[][i]} instead, but beware of differences between matrices and dataframes}
#!      \item{\code{x[cbind()]}}{matrix index subsetting is not defined, use \code{x[][cbind(i, j)]} instead}
#!      \item{\code{ref=TRUE}}{parameter \code{ref} needs to be used sensibly to exploit the advantages of refdata objects}
#!   }
#! }
#! \value{
#!   an object of class refdata (appended to class attributes of data), which is an empty list with two attributes
#!   \item{dat}{the environment where the data x and its dimension dim is stored}
#!   \item{ind}{the environment where the indexes i, j and the effective subset size ni, nj is stored}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{Extract}},  \code{\link[base]{matrix}},  \code{\link[base]{data.frame}}, \code{\link{optimal.index}}, \code{\link{ref}} }
#!
#! \examples{
#!
#!   ## Simple usage Example
#!   x <- cbind(1:5, 5:1)       # take a matrix or data frame
#!   rx <- refdata(x)           # wrap it into an refdata object
#!   rx                         # see the autoprinting
#!   rm(x)                      # delete original to save memory
#!   rx[]                       # extract all data
#!   rx[-1, ]                   # extract part of data
#!   rx2 <- rx[-1, , ref=TRUE]  # create refdata object referencing part of data 
#!                              # (only index, no data is duplicated)
#!   rx2                        # compare autoprinting
#!   rx2[]                      # extract 'all' data
#!   rx2[-1, ]                  # extract part of (part of) data
#!   cat("for more examples look the help pages\n")
#!
#!  \dontrun{
#!   # Memory saving demos
#!   square.matrix.size <- 1000
#!   recursion.depth.limit <- 10
#!   non.referenced.matrix <- matrix(1:(square.matrix.size*square.matrix.size)
#!   , nrow=square.matrix.size, ncol=square.matrix.size)
#!   rownames(non.referenced.matrix) <- paste("a", seq(length=square.matrix.size), sep="")
#!   colnames(non.referenced.matrix) <- paste("b", seq(length=square.matrix.size), sep="")
#!   referenced.matrix <- refdata(non.referenced.matrix)
#!   recurse.nonref <- function(m, depth.limit=10){
#!     x <- m[1,1]   # need read access here to create local copy
#!     gc()
#!     cat("depth.limit=", depth.limit, "  memory.size=", memsize.wrapper(), "\n", sep="")
#!     if (depth.limit)
#!       Recall(m[-1, -1, drop=FALSE], depth.limit=depth.limit-1)
#!     invisible()
#!   }
#!   recurse.ref <- function(m, depth.limit=10){
#!     x <- m[1,1]   # read access, otherwise nothing happens
#!     gc()
#!     cat("depth.limit=", depth.limit, "  memory.size=",  memsize.wrapper(), "\n", sep="")
#!     if (depth.limit)
#!       Recall(m[-1, -1, ref=TRUE], depth.limit=depth.limit-1)
#!     invisible()
#!   }
#!   gc()
#!   memsize.wrapper()
#!   recurse.ref(referenced.matrix, recursion.depth.limit)
#!   gc()
#!    memsize.wrapper()
#!   recurse.nonref(non.referenced.matrix, recursion.depth.limit)
#!   gc()
#!    memsize.wrapper()
#!   rm(recurse.nonref, recurse.ref, non.referenced.matrix
#!   , referenced.matrix, square.matrix.size, recursion.depth.limit)
#!   }
#!   cat("for even more examples look at regression.test.refdata()\n")
#!   regression.test.refdata()  # testing correctness of refdata functionality
#! }
#! \keyword{ programming }
#! \keyword{ manip }

refdata <- function(x){
  d <- dim(x)
  stopifnot(length(d)==2)
  dat <- new.env()
  ind <- new.env()
  assign("x", x, dat)
  assign("dim", d, dat)
  assign("ni", d[1], ind)
  assign("nj", d[2], ind)
  assign("i", NULL, ind)
  assign("j", NULL, ind)
  ref <- list()
  attributes(ref) <- list(dat=dat, ind=ind, class=c("refdata", class(x)))
  ref
}

derefdata <- function(x){
  stopifnot(inherits(x, "refdata"))
  get("x", envir=attr(x,'dat'))
}

"derefdata<-" <- function(x, value){
  stopifnot(inherits(x, "refdata"))
  dat <- attr(x,'dat')
  stopifnot(identical(dim(value), dim(get("x", envir=dat))))
  assign("x", value, envir=dat)
  x
}


"[.refdata" <- function(x, i=NULL, j=NULL, drop=FALSE, ref=FALSE){
  if (!is.null(dim(i)))
    stop("x[cbind(i, j)] matrix subsetting undefined for refdata objects, you can use x[][cbind(i,j)] instead")
  if ( xor(missing(i), missing(j)) && ( nargs() + missing(x) + missing(drop) + missing(ref)) == 4 )
    stop("x[i] single index subsetting undefined for refdata objects, you can use x[][i] instead")
  dat <- attr(x, "dat")
  x2 <- get("x", dat)
  d <- get("dim", dat)
  if (ref){
    new.ind <- new.env()
    ind <- attr(x, "ind")
    if (is.null(i)){
      assign("i", get("i", ind), new.ind)
      assign("ni", get("ni", ind), new.ind)
    }else{
      i <- optimal.index(i, d[1], rownames(x2), i.previous=get("i", ind))
      assign("i", i, new.ind)
      assign("ni", attr(i, "ni"), new.ind)
    }
    if (is.null(j)){
      assign("j", get("j", ind), new.ind)
      assign("nj", get("nj", ind), new.ind)
    }else{
      j <- optimal.index(j, d[2], colnames(x2), i.previous=get("j", ind))
      assign("j", j, new.ind)
      assign("nj", attr(j, "ni"), new.ind)
    }
    attr(x, "ind") <- new.ind
    x
  }else{
    ind <- attr(x, "ind")
    temp <- get("i", ind)
    if (!is.null(temp)){
      if (is.null(i)){
        i <- temp
      }else{
        i <- optimal.index(i, d[1], rownames(x2), i.previous=temp, strict=FALSE)
      }
    }
    temp <- get("j", ind)
    if (!is.null(temp)){
      if (is.null(j)){
        j <- temp
      }else{
        j <- optimal.index(j, d[2], colnames(x2), i.previous=temp, strict=FALSE)
      }
    }
    if (is.null(i)){
      if (is.null(j)){
        x2
      }else{
        x2[, j, drop=drop]
      }
    }else{
      if (is.null(j)){
        x2[i, , drop=drop]
      }else{
        x2[i, j, drop=drop]
      }
    }
  }
}


"[<-.refdata" <- function(x, i=NULL, j=NULL, ref=FALSE, value){
  if (!is.null(dim(i)))
    stop("x[cbind(i, j)] matrix subsetting undefined for refdata objects, you can use x[][cbind(i,j)] instead")
  if ( xor(missing(i), missing(j)) && (nargs() + missing(ref) + missing(value) + missing(x))==4 )
    stop("x[i] single index subsetting undefined for refdata objects, you can use x[][i] instead")
  ind <- attr(x, "ind")
  dat <- attr(x, "dat")
  d <- get("dim", dat)
  if (ref){
    x2 <- get("x", dat)
    if (!is.null(i)){
      i <- optimal.index(i, d[1], rownames(x2), i.previous=get("i", ind))
    }else{
      i <- get("i", ind)
    }
    if (!is.null(j)){
      j <- optimal.index(j, d[2], colnames(x2), i.previous=get("j", ind))
    }else{
      j <- get("j", ind)
    }
    if (is.null(i)){
      if (is.null(j)){
        eval(substitute(x[] <- value, list(value=value)), dat)
      }else{
        eval(substitute(x[, j] <- value, list(j=j, value=value)), dat)
      }
    }else{
      if (is.null(j)){
        eval(substitute(x[i, ] <- value, list(i=i, value=value)), dat)
      }else{
        eval(substitute(x[i, j] <- value, list(i=i, j=j, value=value)), dat)
      }
    }
    x
  }else{
    ii <- get("i", ind)
    jj <- get("j", ind)
    if (is.null(ii)){
      if (is.null(jj)){
        x <- get("x", dat)
      }else{
        x <- get("x", dat)[, jj, drop=FALSE]
      }
    }else{
      if (is.null(jj)){
        x <- get("x", dat)[ii, , drop=FALSE]
      }else{
        x <- get("x", dat)[ii, jj, drop=FALSE]
      }
    }
    if (is.null(i)){
      if (is.null(j)){
        x[] <- value
      }else{
        x[, j] <- value
      }
    }else{
      if (is.null(j)){
        x[i, ] <- value
      }else{
        x[i, j] <- value
      }
    }
    refdata(x)
  }
}



"$.refdata" <- function(x, j, drop=TRUE, ref=FALSE)
  # xx FIXME TODO lazy implementation as special case of [.refdata, can be performance tuned
{
  x[, j, drop=drop, ref=ref]
}

"$<-.refdata" <- function(x, j, ref=FALSE, value)
  # xx FIXME TODO lazy implementation as special case of [<-.refdata, can be performance tuned
{
  x[, j, ref=ref] <- value
}


"[[.refdata" <- function(x, j, drop=TRUE, ref=FALSE)
  # xx FIXME TODO lazy implementation as special case of [.refdata, can be performance tuned
{
  x[, j, drop=drop, ref=ref]
}

"[[<-.refdata" <- function(x, j, ref=FALSE, value)
  # xx FIXME TODO lazy implementation as special case of [<-.refdata, can be performance tuned
{
  x[, j, ref=ref] <- value
}


row.names.refdata <- function(x){
  row.names(x[])
}

names.refdata <- function(x){
  names(x[])
}

dim.refdata <- function(x){
  ind <- attr(x, "ind")
  c(get("ni", ind), get("nj", ind))
}

dimnames.refdata <- function(x){
  dimnames(x[])
}


print.refdata <- function(x, ...){
  dim.dat <- dim(derefdata(x))
  dim.ind <- dim(x)
  cat("refdata (", if (inherits(x, "data.frame")) "data.frame" else if (is.matrix(x[])) "matrix" else "unknown embedded", ") with [", paste(dim.ind, collapse=",") ,"] of [", paste(dim.dat, collapse=","), "]\n", sep="")
  cat("use  x[]  to get the complete actual subset\n")
  cat("use  x[...]  for standard extraction\n")
  cat("use  x[..., ref=TRUE]  to get a newly indexed refdata object\n")
  cat("use  x[...] <- value  to overwrite x with a refdata object containing a new env containing a modified dataset\n")
  cat("use  x[..., ref=TRUE] <- value  to modify the original dataset\n")
}


#! \name{regression.test.refdata}
#! \alias{regression.test.refdata}
#! \title{ regression test for refdata }
#! \description{
#!   This function checks a series of use cases.
#! }
#! \usage{
#! regression.test.refdata()
#! }
#! \details{
#!   raises a stop error if a problem is detected
#! }
#! \value{
#!   TRUE if successful
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{refdata}}, \code{\link[utils]{example}} }
#!
#! \examples{
#!   regression.test.refdata()
#! }
#! \keyword{ internal }


regression.test.refdata <- function(){

  for (i in c("matrix", "data.frame")){

    # example data
    x <- matrix(1:9, ncol=3)
    dimnames(x) <- list(a=paste("a", seq(length=nrow(x)), sep=""), b=paste("b", seq(length=nrow(x)), sep=""))
    if (i=="data.frame"){
      x <- as.data.frame(x)
      cat("testing refdata with data.frame\n")
    }else{
      cat("testing refdata with matrix\n")
    }

    # currently no checks on $.refdata as this is identical to [.refdata

    # check row reduction
    rx3 <- refdata(x)
    stopifnot(identical(rx3[-1, 1:2], x[-1, 1:2]))
    stopifnot(identical(dim(rx3), dim(x)))
    stopifnot(identical(dim(derefdata(rx3)), dim(x)))
    stopifnot(identical(dimnames(rx3), dimnames(x)))
    stopifnot(identical(dimnames(derefdata(rx3)), dimnames(x)))
    rx2 <- rx3[-1, , ref=TRUE]
    stopifnot(identical(rx2[, 1:2], rx3[-1, 1:2]))
    stopifnot(identical(dim(rx2), dim(x[-1, , drop=FALSE])))
    stopifnot(identical(dim(derefdata(rx2)), dim(x)))
    stopifnot(identical(dimnames(rx2), dimnames(x[-1, , drop=FALSE])))
    stopifnot(identical(dimnames(derefdata(rx2)), dimnames(x)))
    rx1 <- rx2[-1, , ref=TRUE]
    stopifnot(identical(rx1[, 1:2], rx2[-1, 1:2]))
    stopifnot(identical(dim(rx1), dim(x[-1, , drop=FALSE][-1, , drop=FALSE])))
    stopifnot(identical(dim(derefdata(rx1)), dim(x)))
    stopifnot(identical(dimnames(rx1), dimnames(x[-1, , drop=FALSE][-1, , drop=FALSE])))
    stopifnot(identical(dimnames(derefdata(rx1)), dimnames(x)))

    # check col reduction
    rx3 <- refdata(x)
    stopifnot(identical(rx3[1:2, -1], x[1:2, -1]))
    rx2 <- rx3[, -1, ref=TRUE]
    stopifnot(identical(rx2[1:2, ], rx3[1:2, -1]))
    stopifnot(identical(dim(rx2), dim(x[, -1, drop=FALSE])))
    stopifnot(identical(dim(derefdata(rx2)), dim(x)))
    stopifnot(identical(dimnames(rx2), dimnames(x[, -1, drop=FALSE])))
    stopifnot(identical(dimnames(derefdata(rx2)), dimnames(x)))
    rx1 <- rx2[, -1, ref=TRUE]
    stopifnot(identical(rx1[1:2, ], rx2[1:2, -1]))
    stopifnot(identical(dim(rx1), dim(x[, -1, drop=FALSE][, -1, drop=FALSE])))
    stopifnot(identical(dim(derefdata(rx1)), dim(x)))
    stopifnot(identical(dimnames(rx1), dimnames(x[, -1, drop=FALSE][, -1, drop=FALSE])))
    stopifnot(identical(dimnames(derefdata(rx1)), dimnames(x)))

    # check row+col reduction
    rx3 <- refdata(x)
    stopifnot(identical(rx3[], x))
    rx2 <- rx3[-1, -1, ref=TRUE]
    stopifnot(identical(rx2[], rx3[-1, -1]))
    rx1 <- rx2[-1, -1, ref=TRUE]
    stopifnot(identical(dim(rx2), dim(x[-1, -1, drop=FALSE])))
    stopifnot(identical(dim(derefdata(rx2)), dim(x)))
    stopifnot(identical(dimnames(rx2), dimnames(x[-1, -1, drop=FALSE])))
    stopifnot(identical(dimnames(derefdata(rx2)), dimnames(x)))
    stopifnot(identical(rx1[], rx2[-1, -1]))
    rx0 <- rx1[-1, -1, ref=TRUE]
    stopifnot(identical(rx0[], rx1[-1, -1]))
    stopifnot(identical(dim(rx1), dim(x[-1, -1, drop=FALSE][-1, -1, drop=FALSE])))
    stopifnot(identical(dim(derefdata(rx1)), dim(x)))
    stopifnot(identical(dimnames(rx1), dimnames(x[-1, -1, drop=FALSE][-1, -1, drop=FALSE])))
    stopifnot(identical(dimnames(derefdata(rx1)), dimnames(x)))

    # check dim dropping
    rx3 <- refdata(x)
    rx1 <- rx3[1, , ref=TRUE]
    stopifnot(identical(rx1[], rx3[1,]))
    stopifnot(identical(rx1[drop=FALSE], rx3[1, , drop=FALSE]))
    stopifnot(identical(rx1[drop=TRUE], rx3[1, , drop=TRUE]))
    stopifnot(identical(rx1[drop=FALSE], x[1, , drop=FALSE]))
    stopifnot(identical(rx1[drop=TRUE], x[1, , drop=TRUE]))
    rx1 <- rx3[, 1, ref=TRUE]
    stopifnot(identical(rx1[], rx3[, 1]))
    stopifnot(identical(rx1[drop=FALSE], rx3[, 1, drop=FALSE]))
    stopifnot(identical(rx1[drop=TRUE], rx3[, 1, drop=TRUE]))
    stopifnot(identical(rx1[drop=FALSE], x[, 1, drop=FALSE]))
    stopifnot(identical(rx1[drop=TRUE], x[, 1, drop=TRUE]))
    rx1 <- rx3[1, 1, ref=TRUE]
    stopifnot(identical(rx1[], rx3[1, 1]))
    stopifnot(identical(rx1[drop=FALSE], rx3[1, 1, drop=FALSE]))
    stopifnot(identical(rx1[drop=TRUE], rx3[1, 1, drop=TRUE]))
    stopifnot(identical(rx1[drop=FALSE], x[1, 1, drop=FALSE]))
    stopifnot(identical(rx1[drop=TRUE], x[1, 1, drop=TRUE]))

    #check assignments
    rx3 <- refdata(x)
    rx2 <- rx3[-1, -1, ref=TRUE]
    rx2b <- rx3[-1, -1, ref=TRUE]
    rx2[] <- x[-1, -1]-1
    stopifnot(identical(rx3[], x))
    stopifnot(identical(rx2[], x[-1, -1] - 1))
    stopifnot(identical(rx2b[], x[-1, -1]))
    rx2 <- rx3[-1, -1, ref=TRUE]
    rx2b <- rx3[-1, -1, ref=TRUE]
    rx2[-1, -1] <- x[-1, -1][-1, -1] - 1
    stopifnot(identical(rx3[], x))
    stopifnot(identical(rx2[-1, -1], x[-1, -1][-1, -1, drop=FALSE] - 1))
    stopifnot(identical(rx2b[-1, -1], x[-1, -1][-1, -1, drop=FALSE]))

    rx2 <- rx3[-1, -1, ref=TRUE]
    rx2b <- rx3[-1, -1, ref=TRUE]
    rx2[, , ref=TRUE] <- x[-1, -1] - 1
    y <- x
    y[-1,-1] <- x[-1, -1] - 1
    stopifnot(identical(rx3[], y))
    stopifnot(identical(rx2[], rx3[-1, -1]))
    stopifnot(identical(rx2b[], rx3[-1, -1]))

  }

}

