# ff dataframe objects
# (c) 2009 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2008-12-31
# Last changed: 2008-12-31

# source("d:/mwp/eanalysis/ff/R/ffdf.R")


if (FALSE){
  library(ff)
  n <- 160
  x <- ff(as.factor(letters[1:16]), dim=c(n,4), dimorder=2:1)
  y <- ff(1:(n*4), dim=c(n,4), dimorder=2:1)
  colnames(y) <- letters[1:4]
  z <- ff(as.factor(letters[1:16]), length=n)
  z2 <- ff(1:n)
  dim(z2) <- c(n,1)
  dimorder <- 2:1
  z3 <- ff(1:n, dim=c(n,1), dimorder=2:1)
  colnames(z3) <- "z3col"
  #rnam <- fffc("", maxwidth=nchar(as.character.hexmode(n)), length=n)
  rnam <- character(n)
  ffvecapply(rnam[i1:i2] <- as.character.hexmode(i1:i2), N=n, VMODE="double", VERBOSE=TRUE)

  #d <- ffdf(x, y, z, z2, z3)
  #dj <- ffdf(x, y, z, z2, z3, ff_join=list(c(1,3), c(2,4,5)), ff_args=list(pattern="dj_"))
  #di <- ffdf(I(x), y, z, z2, z3)
  #dji <- ffdf(I(x), y, z, z2, z3, ff_join=list(c(1,3), c(2,4,5)), ff_args=list(pattern="dj_"))

  d   <- ffdf(x, y, z, z2, z3, row.names=rnam)
  dj  <- ffdf(x, y, z, z2, z3, ff_join=list(c(1,3), c(2,4,5)), row.names=rnam, ff_args=list(pattern="dj_"))
  ds  <- ffdf(x, y, z, z2, z3, ff_split=c(1), row.names=rnam, ff_args=list(pattern="dj_"))
  di  <- ffdf(I(x), y, z, z2, z3, row.names=rnam)
  dji <- ffdf(I(x), y, z, z2, z3, ff_join=list(c(1,3), c(2,4,5)), row.names=rnam, ff_args=list(pattern="dj_"))
  dsi <- ffdf(I(x), y, z, z2, z3, ff_split=c(2), row.names=rnam, ff_args=list(pattern="dj_"))

  d2   <- clone(d  )
  dj2  <- clone(dj )
  ds2  <- clone(ds )
  di2  <- clone(di )
  dji2 <- clone(dji)
  dsi2 <- clone(dsi)

  dj2[1:10,] <- dj[1:10,]

  update(d, d2)
  update(dj, dj2)
  update(ds, ds2)


  i <- quote(1:10)
  d[i,]
  d[i,] <- d[i,]
  d[i,]

  dj[i,]
  dj[i,] <- dj[i,]
  dj[i,]

  ds[i,]
  ds[i,] <- ds[i,]
  ds[i,]

  di[i,]
  di[i,] <- di[i,]
  di[i,]

  dji[i,]
  dji[i,] <- dji[i,]
  dji[i,]

  dsi[i,]
  dsi[i,] <- dsi[i,]
  dsi[i,]


  Rprof(filename = "d:/tmp/d.out", interval = 0.01)
  system.time(x <- d[ii,])
  Rprof(NULL)
  summaryRprof(filename = "d:/tmp/d.out")


  # physical dimorder tests

  library(ff)
  options("fftempdir"="d:/tmp")
  n <- 80000000
  x0 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x1 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x2 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x3 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x4 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x5 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x6 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x7 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x8 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  x9 <- ff(1:n, vmode="integer", caching="mmeachflush", finalizer="close")
  save.image("d:/tmp/ffdf.RData")

  d <- ffdf(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)

  ffrowapply( {
    i <- hi(i1,i2, maxindex=nrow(d))
    #fails: d$x0[i] <- rowSums(d[i, -1, drop=FALSE])
    d[i,1] <- rowSums(d[i, -1, drop=FALSE])
  }
  , X=d
  , BATCHSIZE=1000000
  , VERBOSE=TRUE
  )


  d0 <- colSums(ffrowapply( {
    i <- hi(i1,i2, maxindex=nrow(d))
    colSums(d[i, , drop=FALSE])
  }
  , X=d
  , BATCHSIZE=1000000
  , RETURN = TRUE
  , CFUN = "rbind"
  , VERBOSE=TRUE
  ))
  for (j in 1:10)
  print(ffvecapply( {
    i <- hi(i1,i2, maxindex=nrow(d))
    sum(as.double(d[[j]][i]))
  }
  , X=d[[j]]
  , BATCHSIZE=1000000
  , RETURN = TRUE
  , CFUN = "sum"
  , VERBOSE=TRUE
  ))



  #dj <- ffdf(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, join=list(1:10), ffargs=list(caching="mmeachflush"))
  x <- ff(dim=c(n, 10), vmode="integer", caching="mmeachflush")
  colnames(x) <- paste("x", 0:9, sep="")
  dj <- ffdf(x)

  #emulate partitioning
  y <- vector("list", 10)
  for (h in 1:10){
    message(h, "\n")
    y[[h]] <- ff(1:n, dim=c(n/10,10), dimorder=2:1, vmode="integer", caching="mmeachflush")
  }



  ns <- 10

  i <- sample(n, ns)
  p <- .subset2(d, "physical")
  z <- vector("list", 10)
  system.time({
  for (h in 1:10){
    z[[h]] <- p[[h]][i]
  }
  z <- do.call("cbind", z)
  })
  p <- .subset2(dj, "physical")
  system.time(
    z <- p[[1]][i,]
  )

  i <- i - 1L
  n10 <- n/10L
  z <- vector("list", 10)
  system.time({
    ii <- split(as.integer(i%%(n10))+1L, (i%/%(n10))+1L)
    iin <- as.integer(names(ii))
    for (h in seq_len(length(ii))){
      j <- ii[[h]]
      z[[h]] <- y[[iin[h]]][j,]
    }
    z <- do.call("rbind", z)
  })



  i <- sample(n-ns+1L, 1)
  j <- i+ns-1L
  p <- .subset2(d, "physical")
  system.time(
  for (h in 1:10){
    p[[h]][i:j]
  })
  p <- .subset2(dj, "physical")
  system.time(
    p[[1]][i:j,]
  )
  i <- sample(n/10-ns+1L, 1)
  j <- i+ns-1L
  p <- y[[sample(10, 1)]]
  system.time(
    p[i:j,]
  )


  i <- sample(n, ns)
  system.time(d[i,])
  system.time(dj[i,])

  i <- sample(n, ns)
  j <- sample(10, 5)
  system.time(d[i,j])
  system.time(dj[i,j])

  i <- sample(n, ns)
  system.time(d[i,1])
  system.time(dj[i,1])


  i <- sample(n-ns+1L, 1)
  j <- i+ns-1L
  system.time(for(h in 1:10)d[i:j,])
  system.time(for(h in 1:10)dj[i:j,])

  i <- sample(n-ns+1L, 1)
  j <- i+ns-1L
  system.time(for(h in 1:10)d[i:j,1:5])
  system.time(for(h in 1:10)dj[i:j,1:5])

  i <- sample(n-ns+1L, 1)
  j <- i+ns-1L
  system.time(for(h in 1:10)d[i:j,1])
  system.time(for(h in 1:10)dj[i:j,1])


  Rprof(filename = "d:/tmp/d.out", interval = 0.01)
  x <- d[i,]
  Rprof(NULL)
  rm(x)
  Rprof(filename = "d:/tmp/dj.out", interval = 0.01)
  x <- dj[i,]
  Rprof(NULL)
  rm(x)

  summaryRprof(filename = "d:/tmp/d.out")
  summaryRprof(filename = "d:/tmp/dj.out")
}


#! \name{Forbidden_ffdf}
#! \Rdversion{1.1}
#! \alias{physical<-.ffdf}
#! \alias{virtual<-.ffdf}
#! \alias{length<-.ffdf}
#! \alias{vmode<-.ffdf}
#! \alias{vw<-.ffdf}
#! \alias{vw.ffdf}
#! \title{
#!   Forbidden ffdf functions
#! }
#! \description{
#!   Methods implemented just to prevent using them (because something inppropriate could be find by inheritance)
#! }
#! \usage{
#! \method{physical}{ffdf}(x) <- value
#! \method{virtual}{ffdf}(x) <- value
#! \method{length}{ffdf}(x) <- value
#! \method{vmode}{ffdf}(x, ...) <- value
#! \method{vw}{ffdf}(x, ...) <- value
#! \method{vw}{ffdf}(x, ...)
#! }
#! \arguments{
#!   \item{x}{an ffdf object}
#! }
#! \keyword{ internal }


#! \name{Internal_ffdf}
#! \Rdversion{1.1}
#! \alias{get_nvw}
#! \title{
#!   Internal ffdf functions
#! }
#! \description{
#!   Internal ffdf functions
#! }
#! \usage{
#! get_nvw(x)
#! }
#! \arguments{
#!   \item{x}{an object}
#! }
#! \details{
#! 'get_nvw' gives the length and vw of an object (or its first dimension, if it is not a vector)
#! }
#! \keyword{ internal }

get_nvw <- function(x){
  d <- dim(x)
  if (is.null(d)){
    n <- length(x)
    pvw <- vw(x)
  }else{
    n <- d[[1]]
    pvw <- vw(x)[,1]
  }
  list(n=n, vw=pvw)
}


#! \name{ffdf}
#! \Rdversion{1.1}
#! \alias{ffdf}
#! \title{
#! ff class for data.frames
#! }
#! \description{
#! Function 'ffdf' creates ff data.frames stored on disk very similar to 'data.frame'
#! }
#! \usage{
#! ffdf(...
#! , row.names = NULL
#! , ff_split = NULL
#! , ff_join = NULL
#! , ff_args = NULL
#! , update = TRUE
#! , BATCHSIZE = .Machine$integer.max
#! , BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{\dots}{
#! \code{\link{ff}} vectors or matrices (optionally wrapped in \code{I()} that shall be bound together to an ffdf object
#! }
#!   \item{row.names}{
#! A \code{\link{character}} vector. Not recommended for large objects with many rows.
#! }
#!   \item{ff_split}{
#! A vector of character names or integer positions identifying input components to physically split into single ff_vectors.
#! If vector elements have names, these are used as root name for the new ff files.
#! }
#!   \item{ff_join}{
#! A list of vectors with character names or integer positions identifying input components to physically join in the same ff matrix.
#! If list elements have names, these are used to name the new ff files.
#! }
#!   \item{update}{
#! By default (TRUE) new ff files are updated with content of input ff objects. Setting to FALSE prevents this update.
#! }
#!   \item{ff_args}{
#! a list with further arguments passed to \code{\link{ff}} in case that new ff objects are created via 'ff_split' or 'ff_join'
#! }
#! \item{BATCHSIZE}{
#!   passed to \code{\link{update.ff}}
#! }
#!   \item{BATCHBYTES}{
#!   passed to \code{\link{update.ff}}
#! }
#!   \item{VERBOSE}{
#!   passed to \code{\link{update.ff}}
#! }
#! }
#! \details{
#! By default, creating an 'ffdf' object will NOT create new ff files, instead existing files are referenced.
#! This differs from \code{\link[base]{data.frame}}, which always creates copies of the input objects,
#! most notably in \code{data.frame(matrix())}, where an input matrix is converted to single columns.
#! ffdf by contrast, will store an input matrix physically as the same matrix and virtually map it to columns.
#! Physically copying a large ff matrix to single ff vectors can be expensive.
#! More generally, ffdf objects have a \code{\link[=physical.ffdf]{physical}} and a \code{\link[=virtual.ffdf]{virtual}} component,
#! which allows very flexible dataframe designs: a physically stored matrix can be virtually mapped to single columns,
#! a couple of physically stored vectors can be virtually mapped to a single matrix.
#! The means to configure these are \code{\link[base]{I}} for the virtual representation and the 'ff_split' and 'ff_join'
#! arguments for the physical representation. An ff matrix wrapped into 'I()' will return the input matrix as a single object,
#! using 'ff_split' will store this matrix as single vectors - and thus create new ff files.
#! 'ff_join' will copy a couple of input vectors into a unified new ff matrix with \code{dimorder=c(2,1)},
#! but virtually they will remain single columns. The returned ffdf object has also a \code{\link[=dimorder.ffdf]{dimorder}} attribute,
#! which indicates whether the ffdf object contains a matrix with non-standard dimorder \code{c(2,1)}, see \code{\link{dimorderStandard}}. \cr
#! Currently, \code{\link[=vw]{virtual windows}} are not supported for ffdf.
#! }
#! \section{Methods}{
#!  The following methods and functions are available for ffdf objects:
#!  \tabular{lrll}{
#!   \emph{ Type} \tab  \emph{ Name }  \tab \emph{ Assign }  \tab \emph{Comment}  \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Basic functions}  \cr
#!   function \tab  \code{\link{ffdf}}                       \tab \emph{ }  \tab constructor for ffdf objects \cr
#!   generic  \tab  \code{\link[=update.ffdf]{update}}       \tab \emph{ }  \tab updates one ffdf object with the content of another \cr
#!   generic  \tab  \code{\link[=clone.ffdf]{clone}}         \tab \emph{ }  \tab clones an ffdf object \cr
#!   method   \tab  \code{\link[=print.ffdf]{print}}         \tab \emph{ }  \tab print ffdf \cr
#!   method   \tab  \code{\link[=str.ffdf]{str}}             \tab \emph{ }  \tab ffdf object structure \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Class test and coercion}  \cr
#!   function \tab  \code{\link{is.ffdf}}                    \tab \emph{ }  \tab check if inherits from ff \cr
#!   generic  \tab  \code{\link{as.ffdf}}                    \tab \emph{ }  \tab coerce to ff, if not yet \cr
#!   generic  \tab  \code{\link[=as.data.frame.ffdf]{as.data.frame}}              \tab \emph{ }  \tab coerce to ram data.frame\cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Virtual storage mode} \cr
#!   generic  \tab  \code{\link[=vmode.ffdf]{vmode}}         \tab \code{ }  \tab get virtual modes for all (virtual) columns \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Physical attributes}  \cr
#!   function \tab  \code{\link[=physical.ffdf]{physical}}   \tab \code{ }  \tab get physical attributes \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Virtual attributes} \cr
#!   function \tab  \code{\link[=virtual.ffdf]{virtual}}     \tab \code{ } \tab get virtual attributes \cr
#!   method   \tab  \code{\link[=length.ffdf]{length}}       \tab \code{ }  \tab get length \cr
#!   method   \tab  \code{\link[=dim.ffdf]{dim }}            \tab \code{<-} \tab get dim and set nrow \cr
#!   generic  \tab  \code{\link[=dimorder.ffdf]{dimorder}}   \tab \code{ }  \tab get the dimorder (non-standard if any component is non-standard) \cr
#!   method   \tab  \code{\link[=names.ffdf]{names}}         \tab \code{<-} \tab set and get names \cr
#!   method   \tab  \code{\link[=dimnames.ffdf]{row.names}}  \tab \code{<-} \tab set and get row.names \cr
#!   method   \tab  \code{\link[=dimnames.ffdf]{dimnames}}   \tab \code{<-} \tab set and get dimnames \cr
#!   method   \tab  \code{\link[=filename]{pattern}}         \tab \code{<-} \tab set pattern (rename/move files) \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Access functions}  \cr
#!   method   \tab  \code{\link[=[.ffdf]{[}}                 \tab \emph{<-} \tab set and get data.frame content (\code{[,]}) or get ffdf with less columns (\code{[]}) \cr
#!   method   \tab  \code{\link[=[[.ffdf]{[[}}               \tab \emph{<-} \tab set and get single column ff object \cr
#!   method   \tab  \code{\link[=$.ffdf]{$}}                 \tab \emph{<-} \tab set and get single column ff object \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Opening/Closing/Deleting}                                             \cr
#!   generic  \tab  \code{\link[=is.open.ffdf]{is.open}}     \tab \emph{ }  \tab tri-bool is.open status of the physical ff components \cr
#!   method   \tab  \code{\link[=open.ffdf]{open}}           \tab \emph{ }  \tab open all physical ff objects (is done automatically on access) \cr
#!   method   \tab  \code{\link[=close.ffdf]{close}}         \tab \emph{ }  \tab close all physical ff objects \cr
#!   method   \tab  \code{\link[=delete.ffdf]{delete}}       \tab \emph{ }  \tab deletes all physical ff files \cr
#!   method   \tab  \code{\link[=finalize.ffdf]{finalize}}   \tab \emph{ }  \tab call finalizer \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{processing}                                             \cr
#!   method   \tab  \code{\link[=chunk.ffdf]{chunk}}         \tab \emph{ }  \tab create chunked index \cr
#!   method   \tab  \code{\link[=sortLevels.ffdf]{sortLevels}} \tab \emph{ }  \tab sort and recode levels \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Other}                                                     \cr
#!   }
#! }

#! \value{
#! A list with components
#! \item{physical}{the underlying ff vectors and matrices, to be accessed via \code{\link[=physical.ffdf]{physical}}}
#! \item{virtual}{the virtual features of the ffdf including the virtual-to-physical mapping, to be accessed via \code{\link[=virtual.ffdf]{virtual}}}
#! \item{row.names}{the optional row.names, see argument row.names}
#! and class 'ffdf' (NOTE that ffdf dows not inherit from ff)
#! }
#! \author{
#! Jens Oehlschlägel
#! }
#! \note{
#! Note that in theory, accessing a chunk of rows from a matrix with \code{dimorder=c(2,1)} should be faster than accessing across a bunch of vectors.
#! However, at least under windows, the OS has difficulties filecaching parts from very large files, therefore - until we have partitioning - the recommended physical storage is in single vectors.
#! }
#! \seealso{
#!   \code{\link[base]{data.frame}}, \code{\link{ff}}, for more example see \code{\link[=physical.ffdf]{physical}}
#! }
#! \examples{
#!  m <- matrix(1:12, 3, 4, dimnames=list(c("r1","r2","r3"), c("m1","m2","m3","m4")))
#!  v <- 1:3
#!  ffm <- as.ff(m)
#!  ffv <- as.ff(v)
#!
#!  d <- data.frame(m, v)
#!  ffd <- ffdf(ffm, v=ffv, row.names=row.names(ffm))
#!  all.equal(d, ffd[,])
#!  ffd
#!  physical(ffd)
#!
#!  d <- data.frame(m, v)
#!  ffd <- ffdf(ffm, v=ffv, row.names=row.names(ffm), ff_split=1)
#!  all.equal(d, ffd[,])
#!  ffd
#!  physical(ffd)
#!
#!  d <- data.frame(m, v)
#!  ffd <- ffdf(ffm, v=ffv, row.names=row.names(ffm), ff_join=list(newff=c(1,2)))
#!  all.equal(d, ffd[,])
#!  ffd
#!  physical(ffd)
#!
#!  d <- data.frame(I(m), I(v))
#!  ffd <- ffdf(m=I(ffm), v=I(ffv), row.names=row.names(ffm))
#!  all.equal(d, ffd[,])
#!  ffd
#!  physical(ffd)
#!
#!  rm(ffm,ffv,ffd); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

ffdf <- function(
  ...
, row.names   = NULL
, ff_split    = NULL
, ff_join     = NULL
, ff_args     = NULL
, update      = TRUE
, BATCHSIZE   = .Machine$integer.max
, BATCHBYTES  = getOption("ffbatchbytes")
, VERBOSE     = FALSE
){
  if (!is.null(ff_args$filename))
    stop("use argument 'pattern' (not 'filename') for influencing filenames of ffdf")

  dotstr <- as.character(as.list(substitute(list(...)))[-1])
  ndotstr <- nchar(dotstr)
  dotstr <- ifelse(substr(dotstr, 1, 2)=="I(" & substr(dotstr, ndotstr, ndotstr)==")", substr(dotstr, 3, ndotstr-1), dotstr)
  dotval <- list(...)
  dotlen <- length(dotval)
  if (is.null(names(dotval)))
    dotnam <- character(dotlen)
  else
    dotnam <- names(dotval)
  hasdotnam <- nzchar(dotnam)

  nrows <- NULL
  vnam <- character(0)
  virtual <- NULL
  v <- 0L
  make.virtual <- function(d, p, splitcol=NULL, joining=0){
    f <- dotval[[d]]
    DotIsMatrix <- !is.null(dim(f))
    VirtualIsMatrix <- DotIsMatrix && is.null(splitcol)
    AsIs <- inherits(f, "AsIs")

    if (DotIsMatrix){
      n <- nrow(f)
      fcols <- ncol(f)
      fnam <- colnames(f)
    }else{
      n <- length(f)
      fcols <- 1L
      fnam <- ""
    }
    if (is.null(nrows)){
      nrows <<- n
    }else{
      if (n!=nrows){
        stop("number of rows don't match, recycling not yet implemented")
      }
    }

    if (is.null(splitcol)){

      # don't need this if splitting, only if joining
      if (joining){
        if (is.null(colnames(f)))
          pnam[offset+seq_len(fcols)] <<- ""
        else
          pnam[offset+seq_len(fcols)] <<- colnames(f)
      }

       if (!VirtualIsMatrix || AsIs | inherits(f, "model.matrix")){

         # only one element or treated as one
         v <<- v + 1L
         virtual <<- rbind(virtual, data.frame(
           VirtualVmode        = vmode(f)
         , AsIs                = AsIs
         , VirtualIsMatrix     = VirtualIsMatrix
         , PhysicalIsMatrix    = if (joining>1 || (joining && fcols>1)) TRUE else VirtualIsMatrix
         , PhysicalElementNo   = p
         , PhysicalFirstCol    = offset+1L
         , PhysicalLastCol     = offset+fcols
         , stringsAsFactors    = FALSE
         ))
         vnam[v] <<- if (hasdotnam[d]) dotnam[d] else dotstr[d]
         offset <<- offset + fcols

       }else{

         # matrix input treated as single columns

         if (is.null(fnam)){
           if (hasdotnam[d]){
             fnam <- paste(dotnam[d], seq_len(fcols), sep=".")
           }else{
             fnam <- paste(dotstr[d], seq_len(fcols), sep="")
           }
         }else{
           if (hasdotnam[d])
             fnam <- paste(dotnam[d], fnam, sep=".")
         }
         for (pc in seq_len(fcols)){
           v <<- v + 1L
           virtual <<- rbind(virtual, data.frame(
             VirtualVmode      = vmode(f)
           , AsIs              = AsIs
           , VirtualIsMatrix   = FALSE
           , PhysicalIsMatrix  = TRUE
           , PhysicalElementNo = p
           , PhysicalFirstCol  = offset+1L
           , PhysicalLastCol   = offset+1L
           , stringsAsFactors  = FALSE
           ))
           vnam[v] <<- fnam[pc]
           offset <<- offset + 1L
         }
       }

    }else{
      # case splitting
      if (is.null(fnam)){
        if (hasdotnam[d]){
          fnam <- paste(dotnam[d], splitcol, sep=".")
        }else{
          fnam <- paste(dotstr[d], splitcol, sep="")
        }
      }else{
        if (hasdotnam[d])
          fnam <- paste(dotnam[d], fnam[splitcol], sep=".")
        else
          fnam <- fnam[splitcol]
      }

      v <<- v + 1L
      virtual <<- rbind(virtual, data.frame(
        VirtualVmode        = vmode(f)
      , AsIs                = FALSE
      , VirtualIsMatrix     = FALSE
      , PhysicalIsMatrix    = FALSE
      , PhysicalElementNo   = p
      , PhysicalFirstCol    = offset+1L
      , PhysicalLastCol     = offset+1L
      , stringsAsFactors    = FALSE
      ))
      vnam[v] <<- fnam
      offset <<- offset + 1L
    }
  } # end of make.virtual
  if (is.null(ff_join) && is.null(ff_split)){
    Dimorder <- 1:2  # we reuse the existing physical ff objects, however the virtual layout of the data.frame

    physical <- dotval
    for (p in seq_len(dotlen)){
      offset <- 0L
      f <- dotval[[p]]
      if (!is.ff(f))
        stop("ffdf components must be atomic ff objects")
      make.virtual(d=p, p=p)   # maintaining offset, v, virtual
      if (inherits(f, "AsIs")){
        oldClass(physical[[p]]) <- oldClass(f)[-match("AsIs", oldClass(f))]
      }else{
        dimnames(physical[[p]]) <- NULL
        names(physical[[p]]) <- NULL
      }
    }
    names(physical) <- make.names(ifelse(hasdotnam, dotnam, dotstr), unique=TRUE)

  }else{
    Dimorder <- 2:1  # we create new physical ff objects. Unless update==FALSE, we copy the data to the new objects

    if (is.null(ff_args$pattern))
      ff_args$pattern <- "ffdf"

    ff_args$update <- FALSE
    ff_args$BATCHSIZE <- BATCHSIZE
    ff_args$BATCHBYTES <- BATCHBYTES
    ff_args$VERBOSE <- VERBOSE

    if (!is.null(ff_split)){
      if (is.list(ff_split))
        stop("ff_split must be single vector, not a list")
      splitlen <- length(ff_split)
      if (is.character(ff_split)){
        jsplit <- match(ff_split, names(dotval))
      }else{
        jsplit <- match(ff_split, seq_len(dotlen))
      }
      if (any(is.na(jsplit)))
        stop("could not find ff_split elements", "<", paste(ff_split[is.na(jsplit)], collapse="><"), ">")
      if(any(sapply(jsplit, function(jsi){
        dv <- dotval[[jsi]]
        inherits(dv, "AsIs") || inherits(dv, "model.matrix") || length(dim(dv))!=2
        })
      ))
        stop("can only split matrices that are not I() resp. classes 'AsIs' or 'model.matrix'")
      names(jsplit) <- if (is.null(names(ff_split))) {ifelse(hasdotnam[jsplit], dotnam[jsplit], dotstr[jsplit])} else names(ff_split)
      ff_split <- -jsplit
      rm(jsplit)
    }
    if (!is.null(ff_join)){
      if (!is.list(ff_join))
        stop("ff_join must be a list")
      joinlen <- length(ff_join)
      if (is.null(names(ff_join)))
        names(ff_join) <- rep("join", joinlen)
      for (j in seq_len(joinlen)){
        jjoin <- ff_join[[j]]
        if (is.character(jjoin)){
          jjoin <- match(jjoin, names(dotval))
        }else{
          jjoin <- match(jjoin, seq_len(dotlen))
        }
        if (any(is.na(jjoin)))
          stop("could not find ff_join elements", "<", paste(ff_join[[j]][is.na(jjoin)], collapse="><"), ">")
		ff_join[[j]] <- jjoin
      }
    }
    splitjoin <- c(ff_split, unlist(ff_join))
    if (any(duplicated(splitjoin)))
      stop("some elements duplicatd in ff_split and ff_join")
    ff_rest <- setdiff(seq_len(dotlen), abs(unlist(splitjoin)))
    names(ff_rest) <- ifelse(hasdotnam[ff_rest], dotnam[ff_rest], dotstr[ff_rest])

    # split elements have length=1 and negative position, rest elements are treated like join elements of length 1
    splitjoin <- c(as.list(ff_split), ff_join, as.list(ff_rest))
    splitjoin <- splitjoin[order(sapply(splitjoin, function(sj)abs(sj[[1]])))]
    names(splitjoin) <- make.names(names(splitjoin), unique=TRUE)
    splitjoinlen <- length(splitjoin)

    physical <- vector("list", splitjoinlen)
    p <- 0L
    for (sji in seq_len(splitjoinlen)){
      sjinam <- names(splitjoin)[sji]
      sj <- splitjoin[[sji]]
      nsj <- length(sj)
      if (nsj){
        if (nsj==1 && sj[1]<0){
          # split case
          d <- -sj[1]
          f <- dotval[[d]]
          nf <- ncol(f)
          sjmodes <- vmode(f)
          for (i in seq_len(nf)){
            p <- p + 1L
            offset <- 0L
            make.virtual(d=d, p=p, splitcol=i)  # maintaining offset, v, virtual
            fnam <- paste(sjinam, i, sep=".")
            useffargs <- ff_args
            useffargs[c("vmode", "length", "filename")] <- list(sjmodes, nrows, paste(filename(f), fnam, getOption("ffextension"), sep="."))
            physical[[p]] <- do.call("clone", c(list(f), useffargs))
            names(physical)[p] <- fnam
            if (update){
              vw(f) <- cbind(c(0,nrows,0),c(i-1, 1, nf-i))
              update(physical[[p]], f, BATCHSIZE=BATCHSIZE, BATCHBYTES=BATCHBYTES, VERBOSE=VERBOSE)
              vw(f) <- NULL
            }
          }
        }else{
          # join case including rest
          sjmodes <- character(nsj)
          p <- p + 1L
          offset <- 0L
          pnam <- character(0)
          for (i in seq_len(nsj)){
            d <- sj[i]
            f <- dotval[[d]]
            sjmodes[i] <- vmode(f)
            make.virtual(d=d, p=p, joining=nsj)  # maintaining offset, v, virtual, pnam
          }
          n <- offset  # after returning from all calls to make.virtual we know how many columns in total we got
          if (n==1 && is.null(dim(dotval[[sj[1]]]))){
            # store physically as ff_vector
            f <- dotval[[sj[1]]]
            fnam <- paste(sjinam, "vec.", sep=".")
            useffargs <- ff_args
            useffargs[c("vmode", "length", "pattern")] <- list(sjmodes, nrows, paste(ff_args$pattern, fnam, sep="."))
            physical[[p]] <- do.call("clone", c(list(f), useffargs))
            names(physical)[p] <- fnam
            if (inherits(f, "AsIs")){
              oldClass(f) <- oldClass(f)[-match("AsIs", oldClass(f))]
              names(physical[[p]]) <- names(f)
            }
            if (update){
              update(physical[[p]], f, BATCHSIZE=BATCHSIZE, BATCHBYTES=BATCHBYTES, VERBOSE=VERBOSE)
            }
          }else{
            # store physically as ff_matrix
            offset <- 0L
            anyAsIs <- FALSE
            for (i in seq_len(nsj)){
              f <- dotval[[sj[i]]]
              d <- dim(f)
              if (is.null(d))
                nf <- 1L
              else
                nf <- ncol(f)
              if (inherits(f, "AsIs")){
                oldClass(f) <- oldClass(f)[-match("AsIs", oldClass(f))]
                anyAsIs <- TRUE
              }
              if (i==1){
                fnam <- paste(sjinam, "mat.", sep=".")
                useffargs <- ff_args
                useffargs[c("vmode", "dim","dimorder","pattern")] <- list(names(maxffmode(sjmodes)), c(nrows, n), 2:1, paste(ff_args$pattern, fnam, sep="."))

                physical[[p]] <- do.call("clone", c(list(f), useffargs))
                colnames(physical[[p]]) <- pnam[offset+seq_len(n)]
                names(physical)[p] <- fnam
              }
              if (n!=1 && anyAsIs && is.null(rownames(physical[[p]]))){
                rownames(physical[[p]]) <- rownames(f)
              }
              if (update){
                vw(physical[[p]]) <- cbind(c(0,nrows,0), c(offset, nf, n-offset-nf))
                update(physical[[p]], f, BATCHSIZE=BATCHSIZE, BATCHBYTES=BATCHBYTES, VERBOSE=VERBOSE)
                vw(physical[[p]]) <- NULL
              }
              offset <- offset + nf
            }
          }
        }
      }
    }
  }

  if (!is.null(row.names)){
    #if (! (inherits(row.names, "fffc_vector") || is.character(row.names))  )
    #   stop("row.names(ffdf) must be fffc_vector or character vector")
    if (!is.character(row.names)  )
      stop("row.names(ffdf) must be character vector")
    n <- length(row.names)
    if (is.null(nrows)){
      nrows <- n
    }else{
      if (n!=nrows)
        stop("length(row.names) does not match nrows")
    }
  }

  vnam <- make.names(vnam, unique=TRUE)
  row.names(virtual) <- vnam
  attr(virtual, "Dim") <- c(nrows, v)
  attr(virtual, "Dimorder") <- Dimorder
  ret <- list(
    virtual   = virtual   # list of visible df columns represented as as.integer(c(VirtualIsMatrix, PhysicalIsMatrix, PhysicalElementNo, PhysicalFirstCol, PhysicalLastCol)) where PhysicalFirstCol, PhysicalLastCol are NA if we want all columns = no column selection is required
  , physical  = physical  # list of invisible ff objects
  , row.names = row.names
  )
  oldClass(ret) <- "ffdf"
  ret
}



#! \name{Extract.ffdf}
#! \alias{Extract.ffdf}
#! \alias{[.ffdf}
#! \alias{[<-.ffdf}
#! \alias{[[.ffdf}
#! \alias{[[<-.ffdf}
#! \alias{$.ffdf}
#! \alias{$<-.ffdf}
#! \title{ Reading and writing data.frames (ffdf) }
#! \description{
#!   These are the main methods for reading and writing data from ffdf objects.
#! }
#! \usage{
#! \method{[}{ffdf}(x, i, j, drop = ncols == 1)
#! \method{[}{ffdf}(x, i, j) <- value
#! \method{[[}{ffdf}(x, i, j, exact = TRUE)
#! \method{[[}{ffdf}(x, i, j) <- value
#! \method{$}{ffdf}(x, i)
#! \method{$}{ffdf}(x, i) <- value
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#!   \item{i}{ a row subscript or a matrix subscript or a list subscript }
#!   \item{j}{ a column subscript }
#!   \item{drop}{ logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left. }
#!   \item{value}{ A suitable replacement value: it will be repeated a whole number of times if necessary and it may be coerced: see the Coercion section.  If \code{NULL}, deletes the column if a single column is selected with \code{[[<-} or \code{$<-}. }
#!   \item{exact}{ logical: see \code{\link{[}}, and applies to column names. }
#! }
#! \details{
#!   The subscript methods \code{[}, \code{[[} and \code{$}, behave symmetrical to the assignment functions \code{[<-}, \code{[[<-} and \code{$<-}.
#!   What the former return is the assignment value to the latter.
#!   A notable exception is assigning \code{NULL} in \code{[[<-} and \code{$<-} which removes the \code{\link[=virtual]{virtual}} column from the ffdf (and the \code{\link[=physical]{physical}} component if it is no longer needed by any virtual column).
#!   Creating new columns via \code{[[<-} and \code{$<-} requires giving a name to the new column (character subscripting). \code{[<-} does not allow to create new columns, only to replace existing ones.
#! }
#! \section{Subscript expressions and return values}{
#!   \tabular{rllll}{
#!   \emph{allowed expression}    \tab -- \tab \emph{\code{example}}          \tab -- \tab \emph{\code{returnvalue}}                \cr
#!    row selection  \tab    \tab \code{x[i, ]} \tab    \tab \code{\link{data.frame}} or single row as list if \code{drop=TRUE}, like from data.frame \cr
#!    column selection  \tab    \tab \code{x[ ,i]} \tab    \tab \code{\link{data.frame}} or single column as vector unless \code{drop=TRUE}, like from data.frame  \cr
#!    matrix selection  \tab    \tab \code{x[cbind(i,j)]} \tab    \tab vector of the integer-matrix indexed cells (if the column types are compatible) \cr
#!    virtual selection \tab    \tab \code{x[i]}   \tab    \tab \code{\link{ffdf}}  with the selected columns only \cr
#!    physical selection \tab    \tab \code{x[[i]]} \tab    \tab the selected \code{\link{ff}}            \cr
#!    physical selection \tab    \tab \code{x$i} \tab    \tab the selected \code{\link{ff}}            \cr
#!   }
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ffdf}}, \code{\link[=[.data.frame]{Extract.data.frame}}, \code{\link{Extract.ff}}  }
#! \examples{
#!    d <- data.frame(a=letters, b=rev(letters), c=1:26)
#!    x <- as.ffdf(d)
#!
#!    d[1,]
#!    x[1,]
#!
#!    d[1:2,]
#!    x[1:2,]
#!
#!    d[,1]
#!    x[,1]
#!
#!    d[,1:2]
#!    x[,1:2]
#!
#!    d[cbind(1:2,2:1)]
#!    x[cbind(1:2,2:1)]
#!
#!    d[1]
#!    x[1]
#!
#!    d[[1]]
#!    x[[1]]
#!
#!    d$a
#!    x$a
#!
#!    d$a[1:2]
#!    x$a[1:2]
#!
#!    rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }




# for x[[i,j]]  we delegate to x[i,j]   which is more consistent when element i is a I(matrix)
"[[.ffdf" <- function(x, i, j, exact = TRUE){
  virtual <- .subset2(x,"virtual")
  physical <- .subset2(x,"physical")
  # the dummy-thing here make sure we get the same behaviour as with real lists or data.frames
  r <- row.names(virtual)
  names(r) <- r
  dummy <- as.list(r)
  vnam <- .subset2(dummy, i, exact=exact)
  if (is.null(vnam))
    if (missing(j))
      return(NULL)
    else
      stop("subscript out of bounds")
  else{
    if (!missing(j))
      return(x[vnam,j])
    v <- virtual[vnam, , drop=TRUE]
    p <- physical[[v$PhysicalElementNo]]
    if (v$PhysicalIsMatrix){
      n <- v$PhysicalLastCol - v$PhysicalFirstCol + 1L
      d <- dim(p)
      if (d[[2]]!=n){
        vw <- vw(p)
        if (is.null(vw)){
          vw_low <- 0L
          vw_high <- 0L
        }else{
          vw_low <- vw[1,2]
          vw_high <- vw[3,2]
        }
        low <- v$PhysicalFirstCol-1
        vw(p) <- cbind(c(0, d[[1]], 0), c(vw_low + low, n, d[[2]] - low - n + vw_high))
      }
      if (!missing(j))
        p <- p[,j]
    }else{
      if (!missing(j))
        p <- p[j]
    }
  }
  if (v$AsIs)
    oldClass(p) <- c("AsIs", oldClass(p))
  return(p)
}


# for now we here mimic the bad behaviour of R (partial matching is a pain)
"$.ffdf" <- function(x, i){
  x[[i, exact=FALSE]]
}


# for x[[i,j]] <- value   we delegate to x[i,j] <- value   which is more consistent when element i is a I(matrix)
"[[<-.ffdf" <- function(x, i, j, value){
  d <- dim(x)
  virtual <- .subset2(x,"virtual")
  physical <- .subset2(x,"physical")
  r <- row.names(virtual)
  names(r) <- r
  dummy <- as.list(r)
  vnam <- .subset2(dummy, i, exact=TRUE)

  if (is.null(vnam)){
    if (!missing(j))
      stop("subscript out of bounds")
    if (is.null(value))
      return(x)
    if (!is.character(i))
      stop("new ffdf elements need to be addressed by name")
    dummy[i] <- i
    vnam <- i
    d[[2]] <- d[[2]] + 1L
  }else{
    if (!missing(j)){
      if (is.ff(value) || is.list(value))
        stop("when using 2nd index j the assigned value must be a (scalar) vector")
      else{
        x[i,j] <- value
        return(x)
      }
    }
    # remove old slot if necessary
    vpos <- match(vnam, r)
    oldpos <- unique(virtual$PhysicalElementNo[-vpos])
    virtual$PhysicalElementNo <- match(virtual$PhysicalElementNo, oldpos)
    physical <- physical[oldpos]
  }

  if (is.null(value)){
    virtual <- virtual[-vpos, ]
    d[[2]] <- d[[2]] - 1L
  }else{
    # check value (here we know missing(j)==TRUE)
    if (!is.ff(value))
      stop("assigned value must be ff")
    valuedim <- dim(value)
    if (is.null(valuedim)){
      if (length(value) != d[[1]])
        stop("length of ff_vector must match number of ffdf rows")
    }else{
      if (length(valuedim)!=2)
        stop("ff matrices are the only ff_arrays that can be assigned")
      if (valuedim[[1]] != d[[1]])
        stop("number of rows of ff_matrix must match number of ffdf rows")
    }

    # add new slot
    p <- length(physical)+1L
    v <- list(
      VirtualVmode        = vmode(value)
    , AsIs                = inherits(value, "AsIs")
    , VirtualIsMatrix     = !is.null(valuedim)
    , PhysicalIsMatrix    = !is.null(valuedim)
    , PhysicalElementNo   = p
    , PhysicalFirstCol    = 1L
    , PhysicalLastCol     = if (is.null(valuedim)) 1L else valuedim[[2]]
    )
    if (v$AsIs)
      oldClass(value) <- oldClass(value)[-match("AsIs", oldClass(value))]

    virtual[vnam,] <- v
    pnam <- make.names(c(names(physical), vnam), unique=TRUE)
    physical[p] <- list(value)
    names(physical) <- pnam
  }

  attr(virtual, "Dim") <- d
  attr(virtual, "Dimorder") <- if (all(sapply(physical, function(x)dimorderStandard(dimorder(x))))) 1:2 else 2:1

  cl <- oldClass(x)
  oldClass(x) <- NULL
  x$virtual <- virtual
  x$physical <- physical
  oldClass(x) <- cl
  return(x)
}

"$<-.ffdf" <- function(x, i, value){
  x[[i]] <- value
  x
}




"[.ffdf" <- function (x, i, j, drop = ncols == 1)
{
    d <- dim(x)
    mdrop <- missing(drop)
    Narg <- nargs() - (!mdrop)
    if (Narg < 3) {
        if (!mdrop)
            warning("drop argument will be ignored")
        if (missing(i))
            return(x)
        if (is.matrix(i)){
          ii <- split(i[,1], i[,2])
          ret <- lapply(names(ii), function(i){
            x[[as.integer(i)]][ii[[i]]]
          })
          return(unsplit(ret, i[,2]))
        }

      # Now handling one-dimensional index: return a ffdf object (colum selection thereof)

      cl <- oldClass(x)
      oldClass(x) <- NULL

      virtual <- x$virtual[i, , drop=FALSE]
      if (any(is.na(virtual$VirtualIsMatrix)))
        stop("non-existing column selected")

      oldpos <- unique(virtual$PhysicalElementNo)
      virtual$PhysicalElementNo <- match(virtual$PhysicalElementNo, oldpos)
      physical <- x$physical[oldpos]

      ncols <- nrow(virtual)

      attr(virtual, "Dim") <- c(d[[1]], ncols)
      attr(virtual, "Dimorder") <- if (all(sapply(physical, function(x)dimorderStandard(dimorder(x))))) 1:2 else 2:1

      x$virtual <- virtual
      x$physical <- physical
      oldClass(x) <- cl

      return(x)

    }else{

      ncols <- length(x)

      # first handle i=ff subscript
      if (!missing(i) && is.ff(i)){
        #first reduce number of columns if possible
        if (!missing(j)){
          x <- x[j]
        }
        if (drop && ncol(x)==1){
          return(ffindexget(x[[1]], i))
        }else{
          return(ffdfindexget(x, i))
        }
      }

      rownam <- .subset2(x,"row.names")
      physical <- .subset2(x,"physical")
      if (ncols){

        if (missing(j))
          virtual <- .subset2(x,"virtual")
        else{
          virtual <- .subset2(x,"virtual")[j, , drop=FALSE]
          if (any(is.na(virtual$VirtualIsMatrix)))
            stop("non-existing column selected")
        }
        ncols <- nrow(virtual)
        df <- vector("list", ncols)
        pv <- split(seq_len(ncols), virtual$PhysicalElementNo)
        np <- length(pv)
        last_nvw <- NULL
        for (ip in seq_len(np)){
          v <- pv[[ip]]
          vi <- virtual[v,]
          vi1 <- vi[1,,drop=TRUE]
          p <- physical[[vi1$PhysicalElementNo]]
          nvw <- get_nvw(p)
          if (is.null(last_nvw) || !identical(last_nvw, nvw)){
            last_nvw <- nvw
            # note that the hi index stores positions AFTER translating relative index positions (relative to vw) to absolute index positions
            # and thus index may differ for different physical components of the dataframe (if their vw differ)
            # if you pass-in a hi object, this must be suitable for all physical components
            if (missing(i))
              i2 <- hi(from=1, to=nvw$n, maxindex=nvw$n, vw=nvw$vw, pack=FALSE)
            else{
              i2 <- as.hi(i, maxindex=nvw$n, vw=nvw$vw, pack=FALSE, envir=parent.frame(), names=rownam)
            }
          }
          if (ip==1){
            nrows <- length(i2)
          }else{
            if (length(i2)!=nrows)
              stop("number of rows don't match, recycling not yet implemented")
          }

          if (vi1$PhysicalIsMatrix){
            cols <- cumsum(vi$PhysicalLastCol - vi$PhysicalFirstCol + 1L)
            cols <- vecseq(c(1L, cols[-length(cols)] + 1L), cols, concat=FALSE)
            colindex <- vecseq(vi$PhysicalFirstCol, vi$PhysicalLastCol, eval=FALSE)
            pvalue <- p[i2, colindex, drop=FALSE]
            for (iv in seq_len(length(v))){
              elem <- pvalue[, cols[[iv]], drop=!vi$VirtualIsMatrix[iv]]
              if (vi$AsIs[iv]){
                oldClass(elem) <- c("AsIs", oldClass(elem))
              }else{
                dimnames(elem) <- NULL
                names(elem) <- NULL
              }

              df[[v[iv]]] <- elem
            }
          }else{
            pvalue <- p[i2]
            df[v] <- list(pvalue)
          }
        }
        if (ncols==1 && drop){
          df <- df[[1]]
        }else{
          names(df) <- row.names(virtual)
          if (nrows!=1 || !drop){
            class(df) <- "data.frame"
          }
        }
      }else{
        p <- physical[[1]]
        last_nvw <- nvw <- get_nvw(p)
        if (missing(i)){
          nrows <- nvw$n
        }else{
          i2 <- as.hi(i, maxindex=nvw$n, vw=nvw$vw, pack=FALSE, envir=parent.frame(), names=rownam)
          nrows <- length(i2)
        }
        df <- list()
        class(df) <- "data.frame"
      }
      if (is.data.frame(df)){
        if (is.null(rownam)){
          if (missing(i))
            row.names(df) <- seq_len(nrows)
          else{
						if (bit:::anyDuplicated.rlepack(i2$x))
							row.names(df) <- make.unique(as.character(as.which(i2)))
						else
							row.names(df) <- as.which(i2)
					}
        }else{
          if (is.ff(rownam)){
            nvw <- get_nvw(rownam)
            if (!identical(last_nvw, nvw)){
              if (missing(i))
                i2 <- hi(from=1, to=nvw$n, maxindex=nvw$n, vw=nvw$vw, pack=FALSE)
              else{
                i2 <- as.hi(i, maxindex=nvw$n, vw=nvw$vw, pack=FALSE, envir=parent.frame(), names=rownam)
              }
            }
						if (bit:::anyDuplicated.rlepack(i2$x))
							row.names(df) <- make.unique(as.character(rownam[i2]))
						else
							row.names(df) <- rownam[i2]
          }else{
            if (missing(i)){
              row.names(df) <- rownam
            }else if(is.character(i))
              row.names(df) <- make.unique(i)
            else{
							if (bit:::anyDuplicated.rlepack(i2$x))
								row.names(df) <- make.unique(rownam[as.integer(i2)])
							else
								row.names(df) <- rownam[as.integer(i2)]
            }
          }
        }
      }

      return(df)
    }
}


"[<-.ffdf" <- function (x, i, j, value)
{
    d <- dim(x)
    Narg <- nargs()
    if (Narg < 4) {
      if (missing(i))
          i <- 1:d[[2]]
      else if (is.matrix(i)){
        value <- split(rep(value, length.out=nrow(i)), i[,2])
        ii <- split(i[,1], i[,2])
        lapply(seq_along(ii), function(i){
          x[[i]][ii[[i]]] <- value[[i]]
        })
        return(x)
      }

      # Now handling one-dimensional index: assign ffdf object to (colum selection) of ffdf
      n <- length(i)
      if (!n)
        return(x)

      if (!inherits(value, "ffdf"))
        stop("value must be ffdf if only one index used")
        dv <- dim(value)
      if (dv[[1]] != d[[1]])
        stop("number of rows does not match")
      if (dv[[2]] > n)
        stop("too many columns in assign value")
      if (length(i)>d[2])
        stop("too many columns selected for replacement")

      selector <- seq_len(d[2])
      names(selector) <- colnames(x)
      i <- selector[i]
      if (any(is.na(i)))
        stop("non-existing column selected")

      if (any(duplicated(i)))
        stop("duplicated assgignments not allowed")

      cl <- oldClass(x)
      oldClass(x) <- NULL
      oldClass(value) <- NULL

      virtual <- x$virtual
      oldphy <- virtual$PhysicalElementNo[-i]
      uniqueoldpy <- unique(oldphy)
      virtual$PhysicalElementNo[-i] <- match(oldphy, uniqueoldpy)
      physical <- c(x$physical[uniqueoldpy], value$physical)
      value$virtual$PhysicalElementNo <- length(uniqueoldpy) + value$virtual$PhysicalElementNo

      # recycle value
      j <- repfromto(seq_len(dv[[2]]), 1, n)
      virtual[i,] <- value$virtual[j,]

      names(physical) <- make.names(names(physical), unique=TRUE)

      attr(virtual, "Dimorder") <- if (all(sapply(physical, function(x)dimorderStandard(dimorder(x))))) 1:2 else 2:1

      x$virtual <- virtual
      x$physical <- physical
      oldClass(x) <- cl

      return(x)

    }else{

      ncols <- length(x)

      # first handle i=ff subscript
      if (!missing(i) && is.ff(i)){
        #first reduce number of columns if possible
        if (missing(j))
          y <- x
        else
          y <- x[j]
        if (ncol(y)==1 && is.ff(value)){
          ffindexset(y[[1]], i, value)
          return(x)
        }else if(is.ffdf(value)){
          ffdfindexset(y, i, value)
          return(x)
        }else{
          stop("ff/ffdf-iness of value and selected columns don't match")
        }
      }

      rownam <- .subset2(x,"row.names")
      physical <- .subset2(x,"physical")
      if (ncols){

        valuedim <- dim(value)

        if (missing(j))
          virtual <- .subset2(x,"virtual")
        else{
          virtual <- .subset2(x,"virtual")[j, , drop=FALSE]
          if (any(is.na(virtual$VirtualIsMatrix)))
            stop("non-existing column selected")
        }
        ncols <- nrow(virtual)
        #not needed here: valcols <- cumsum(virtual$PhysicalLastCol - virtual$PhysicalFirstCol + 1L)
        #not needed here: valcols <- vecseq(c(1L, valcols[-length(valcols)] + 1L), valcols, concat=FALSE)

        pv <- split(seq_len(ncols), virtual$PhysicalElementNo)
        np <- length(pv)
        last_nvw <- NULL
        for (ip in seq_len(np)){
          v <- pv[[ip]]
          vi <- virtual[v,]
          vi1 <- vi[1,,drop=TRUE]
          p <- physical[[vi1$PhysicalElementNo]]
          nvw <- get_nvw(p)
          if (is.null(last_nvw) || !identical(last_nvw, nvw)){
            last_nvw <- nvw
            # note that the hi index stores positions AFTER translating relative index positions (relative to vw) to absolute index positions
            # and thus index may differ for different physical components of the dataframe (if their vw differ)
            # if you pass-in a hi object, this must bu suitable for all physical components
            if (missing(i))
              i2 <- hi(from=1, to=nvw$n, maxindex=nvw$n, vw=nvw$vw, pack=FALSE)
            else{
              i2 <- as.hi(i, maxindex=nvw$n, vw=nvw$vw, pack=FALSE, envir=parent.frame(), names=rownam)
            }
          }
          if (ip==1){
            nrows <- length(i2)
            if (is.null(valuedim)){
              valuelen <- length(value)
              if (nrows==1){
                if (ncols %% valuelen)
                  stop("ncol(index) not a multiple of length(value)")
                rft <- repfromto(1:valuelen, 1, ncols)
                f <- function(val, ind)val[rft[ind]]
              }else{
                if (nrows %% valuelen)
                  stop("nrow(index) not a multiple of length(value)")
                f <- function(val, ind)val
              }
            }else{
              if (nrows %% valuedim[[1]])
                stop("nrow(index) not a multiple of nrow(value)")
              if (ncols %% valuedim[[2]]){
                stop("ncol(index) not a multiple of ncol(value)")
              }
              rft <- repfromto(1:valuedim[[2]], 1, ncols)
              if (is.matrix(value))
                f <- function(val, ind)val[,rft[ind]]
              else
                f <- function(val, ind)if (length(ind)>1)as.matrix(val[,rft[ind],drop=FALSE]) else val[,rft[ind]]
            }
          }else{
            if (length(i2)!=nrows)
              stop("number of rows in ffdf don't match other number of rows in ffdf, should not happen")
          }
          if (vi1$PhysicalIsMatrix){
            #not needed here: cols <- cumsum(vi$PhysicalLastCol - vi$PhysicalFirstCol + 1L)
            #not needed here: cols <- vecseq(c(1L, cols[-length(cols)] + 1L), cols, concat=FALSE)
            colindex <- vecseq(vi$PhysicalFirstCol, vi$PhysicalLastCol, eval=FALSE)

            # finally the assignment
            p[i2, colindex] <- f(value, v)
          }else{
            p[i2] <- f(value,v)
          }
        }
      }

      return(x)
    }
}


#! \name{clone.ffdf}
#! \Rdversion{1.1}
#! \alias{clone.ffdf}
#! \title{
#! Cloning ffdf objects
#! }
#! \description{
#! clone physically duplicates ffdf objects
#! }
#! \usage{
#! \method{clone}{ffdf}(x, nrow=NULL, ...)
#! }
#! \arguments{
#!   \item{x}{an \code{\link{ffdf}} }
#!   \item{nrow}{ optionally the desired number of rows in the new object. Currently this works only together with \code{initdata=NULL} }
#!   \item{\dots}{ further arguments passed to \code{\link{clone}} (usually not usefull) }
#! }
#! \details{
#!   Creates a deep copy of an ffdf object by cloning all \code{\link[=physical.ffdf]{physical}} components including the \code{\link[=dimnames.ffdf]{row.names}}
#! }
#! \value{
#!   An object of type \code{\link{ffdf}}
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{clone}}, \code{\link{ffdf}}
#! }
#! \examples{
#!   x <- as.ffdf(data.frame(a=1:26, b=letters))
#!
#!   message("Here we change the content of both x and y by reference")
#!   y <- x
#!   x$a[1] <- -1
#!   y$a[1]
#!
#!   message("Here we change the content only of x because y is a deep copy")
#!   y <- clone(x)
#!   x$a[2] <- -2
#!   y$a[2]
#!   rm(x, y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

clone.ffdf <- function(x, nrow=NULL, ...){
  cl <- oldClass(x)
  oldClass(x) <- NULL

  if (is.null(nrow)){

    x$physical <- lapply(x$physical, clone, ...)
    if (is.ff(x))
        x$row.names <- clone(x$row.names, ...)

  }else{

    nrow <- as.integer(nrow)
    x$physical <- lapply(x$physical, function(o, ...){
        d <- dim(o)
        if (is.null(d))
          clone(o, length=nrow, ...)
        else
          clone(o, dim=c(nrow, d[[2]]), ...)
    }, ...)
    if (is.ff(x)){
      if (is.null(nrow))
        x$row.names <- clone(x$row.names, ...)
      else
        x$row.names <- clone(x$row.names, length=nrow, ...)
    }

    attr(x$virtual, "Dim") <- c(nrow, attr(x$virtual, "Dim")[[2]])
  }


  oldClass(x) <- cl
  x
}




update.ffdf <- function(object, from, ...){
  dobject <- dim(object)
  dfrom <- dim(from)
  if (is.null(dfrom) || length(dfrom)!=2 || dobject[[2]]!=dfrom[[2]])
    stop("from must be also ffdf with the same number of columns in update.ffdf(object, from, ...)")
  if (dobject[[1]]%%dfrom[[1]])
    stop("nrow(object) not a multiple of nrow(from) in update.ffdf(object, from, ...)")
  cellbytes <- mean(.rambytes[vmode(object)])
  i1 <- i2 <- 0L # keep R CMD CHECK quiet
  if (dfrom[[1]]<dobject[[1]]){
    ffrowapply(object[i1:i2,] <- from[repfromto(seq_len(dfrom[[1]]), i1, i2),,drop=FALSE], X=object, VBYTES=cellbytes, ...)
  }else
    ffrowapply(object[i1:i2,] <- from[i1:i2,,drop=FALSE], X=object, VBYTES=cellbytes, ...)
  object
}



#! \name{is.ffdf}
#! \alias{is.ffdf}
#! \title{ Test for class ff }
#! \description{
#!   checks if x inherits from class "ffdf"
#! }
#! \usage{
#! is.ffdf(x)
#! }
#! \arguments{
#!   \item{x}{ any object }
#! }
#! \value{
#!   logical scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{inherits}}, \code{\link{as.ffdf}}, \code{\link{is.ff}} }
#! \examples{
#!   is.ffdf(integer())
#! }
#! \keyword{ IO }
#! \keyword{ data }


is.ffdf <- function(x)
  inherits(x, "ffdf")



#! \name{as.ffdf}
#! \Rdversion{1.1}
#! \alias{as.ffdf}
#! \alias{as.ffdf.ff_vector}
#! \alias{as.ffdf.ff_matrix}
#! \alias{as.ffdf.data.frame}
#! \alias{as.data.frame.ffdf}
#! \title{
#! Coercing to ffdf and data.frame
#! }
#! \description{
#!   Functions for coercing to ffdf and data.frame
#! }
#! \usage{
#! as.ffdf(x, ...)
#! \method{as.ffdf}{ff_vector}(x, ...)
#! \method{as.ffdf}{ff_matrix}(x, ...)
#! \method{as.ffdf}{data.frame}(x, vmode=NULL, col_args = list(), ...)
#! \method{as.data.frame}{ffdf}(x, ...)
#! }
#! \arguments{
#!   \item{x}{ the object to be coerced }
#!   \item{vmode}{ optional specification of the \code{\link{vmode}s} of columns of the \code{\link{data.frame}}. Either a character vector of vmodes (named with column names of the data.frame or recycled if not named)
#!                 or a list named with vmodes where each element identifies those columns of the data.frame that should get the vmode encoded in the name of the element }
#!   \item{col_args}{ further arguments; passed to \code{\link{ff}}  }
#!   \item{\dots}{ further arguments; passed to \code{\link{ffdf}} for .ff_vector, .ff_matrix and .data.frame methods, ignored for .ffdf identity method }
#! }
#! \value{
#!   'as.ffdf' returns an object of class \code{\link{ffdf}}, 'as.data.frame' returns an object of class \code{\link{data.frame}}
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{is.ffdf}}, \code{\link{ffdf}}, \code{\link{data.frame}}
#! }
#! \examples{
#!   d <- data.frame(x=1:26, y=letters, z=Sys.time()+1:26)
#!   ffd <- as.ffdf(d)
#!   stopifnot(identical(d, as.data.frame(ffd)))
#!   rm(ffd); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

as.ffdf.ff_matrix <- function(x, ...){
  ffdf(x, row.names=rownames(x), ...)
}
as.ffdf.ff_vector <- function(x, ...){
  ffdf(x, row.names=names(x), ...)
}
as.ffdf.data.frame <- function(
  x
, vmode = NULL
, col_args=list()
, ...
){
  rnam <- attr(x, "row.names")
  if (is.integer(rnam)){
    if (all(rnam==seq_along(rnam)))
      rnam <- NULL
    else
      rnam <- as.character(rnam)
  }
  x <- as.list(x)
  vmodes <- vector("list", length(x))
  if (!is.null(vmode)){
    nam <- names(x)
    if (is.list(vmode)){
      vnam <- names(vmode)
      i <- match(vnam, .vmode[.vimplemented])
      if (any(is.na(i)))
        stop("vmodes not implemented: ", paste("'", vnam[is.na(i)], "'" , collapse=",", sep=""))
      names(vmodes) <- nam
      for (v in vnam)
        vmodes[vmode[[v]]] <- v
    }else{
      vnam <- names(vmode)
      if (is.null(vnam)){
        vmodes <- as.list(rep(vmode, length.out=length(x)))
      }else{
        i <- match(vnam, nam)
        if (any(is.na(i)))
          stop("names not matched for vmode specification: ", paste("'", vnam[is.na(i)], "'", collapse=",", sep=""))
        vmodes[i] <- as.list(vmode)
      }
    }
  }
  if (is.null(col_args$pattern))
    col_args$pattern <- "ffdf"

  l <- list(...)
  if (is.null(l$ff_args$pattern))
    l$ff_args$pattern <- "ffdf"
  ret <- lapply(seq_along(x)
  , function(i, ...){
    xi <- x[[i]]
    AsIs <- inherits(xi, "AsIs")
    if (AsIs){
      oldClass(xi) <- oldClass(xi)[-match("AsIs", oldClass(xi))]
      ret <- do.call("as.ff", c(list(xi, vmode=vmodes[[i]]), col_args))
      oldClass(ret) <- c("AsIs", oldClass(ret))
      ret
    }else{
      do.call("as.ff", c(list(xi, vmode=vmodes[[i]]), col_args))
    }
  }
  , ...
  )
  names(ret) <- names(x)
  do.call("ffdf", c(ret, list(row.names=rnam), l))
}

as.data.frame.ffdf <- function(x, ...)
  x[,]


#! \name{length.ffdf}
#! \Rdversion{1.1}
#! \alias{length.ffdf}
#! \title{
#! Getting length of a ffdf dataframe
#! }
#! \description{
#! Getting "length" (number of columns) of a ffdf dataframe
#! }
#! \usage{
#! \method{length}{ffdf}(x)
#! }
#! \arguments{
#!   \item{x}{an \code{\link{ffdf}} object}
#! }
#! \value{ integer number of columns}
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{dim.ffdf}}, \code{\link{length.ff}}, \code{\link{ffdf}}
#! }
#! \examples{
#!   length(as.ffdf(data.frame(a=1:26, b=letters)))
#!   gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


length.ffdf <- function(x){
  attr(.subset2(x,"virtual"), "Dim")[[2]]  # avoid unclass(x)$virtual
}

"length<-.ffdf" <- function(x, value){
  stop("not allowed to change the number of (virtual) columns of ffdf")
}


#! \name{vmode.ffdf}
#! \Rdversion{1.1}
#! \alias{vmode.ffdf}
#! \title{
#! Virtual storage mode of ffdf
#! }
#! \description{
#! Function vmode returns the virtual storage mode of each ffdf column
#! }
#! \usage{
#! \method{vmode}{ffdf}(x, ...)
#! }
#! \arguments{
#!   \item{x}{\code{\link{ffdf}}}
#!   \item{\dots}{ignored}
#! }
#! \value{
#!   a character vector with one element for each column
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{ \code{\link{vmode}}, \code{\link{ffdf}} }
#! \examples{
#!   vmode(as.ffdf(data.frame(a=as.double(1:26), b=letters)))
#!   gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


vmode.ffdf <- function(x, ...){
  ret <- sapply(.subset2(x,"physical"), vmode)[.subset2(x,"virtual")$PhysicalElementNo]
  names(ret) <- names(x)
  ret
}

"vmode<-.ffdf" <- function(x, ..., value){
  stop("ffdf objects cannot be coerced to a different vmode")
}






#! \name{chunk.ffdf}
#! \Rdversion{1.1}
#! \alias{chunk.ffdf}
#! \alias{chunk.ff_vector}
#! \title{
#!    Chunk ff_vector and ffdf
#! }
#! \description{
#!    Chunking method for ff_vector and ffdf objects (row-wise) automatically considering RAM requirements from recordsize as calculated from \code{\link{sum}(\link{.rambytes}[\link{vmode}])}
#! }
#! \usage{
#! \method{chunk}{ff_vector}(x
#! , RECORDBYTES = .rambytes[vmode(x)], BATCHBYTES = getOption("ffbatchbytes"), \dots)
#! \method{chunk}{ffdf}(x
#! , RECORDBYTES = sum(.rambytes[vmode(x)]), BATCHBYTES = getOption("ffbatchbytes"), \dots)
#! }
#! \arguments{
#!   \item{x}{\code{\link{ff}} or \code{\link{ffdf}}}
#!   \item{RECORDBYTES}{ optional integer scalar representing the bytes needed to process an element of the \code{ff_vector} a single row of the \code{ffdf} }
#!   \item{BATCHBYTES}{ integer scalar limiting the number of bytes to be processed in one chunk, default from \code{getOption("ffbatchbytes")}, see also \code{\link{.rambytes}} }
#!   \item{\dots}{further arguments passed to \code{\link[bit]{chunk}}}
#! }
#! \value{
#!   A list with \code{\link[bit]{ri}} indexes each representing one chunk
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{ \code{\link[bit]{chunk}}, \code{\link{ffdf}} }
#! \examples{
#!   x <- data.frame(x=as.double(1:26), y=factor(letters), z=ordered(LETTERS))
#!   a <- as.ffdf(x)
#!   ceiling(26 / (300 \%/\% sum(.rambytes[vmode(a)])))
#!   chunk(a, BATCHBYTES=300)
#!   ceiling(13 / (100 \%/\% sum(.rambytes[vmode(a)])))
#!   chunk(a, from=1, to = 13, BATCHBYTES=100)
#!   rm(a); gc()
#!
#!   message("dummy example for linear regression with biglm on ffdf")
#!   library(biglm)
#!
#!   message("NOTE that . in formula requires calculating terms manually
#!     because . as a data-dependant term is not allowed in biglm")
#!   form <- Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species
#!
#!   lmfit <- lm(form, data=iris)
#!
#!   firis <- as.ffdf(iris)
#!   for (i in chunk(firis, by=50)){
#!     if (i[1]==1){
#!       message("first chunk is: ", i[[1]],":",i[[2]])
#!       biglmfit <- biglm(form, data=firis[i,,drop=FALSE])
#!     }else{
#!       message("next chunk is: ", i[[1]],":",i[[2]])
#!       biglmfit <- update(biglmfit, firis[i,,drop=FALSE])
#!     }
#!   }
#!
#!   summary(lmfit)
#!   summary(biglmfit)
#!   stopifnot(all.equal(coef(lmfit), coef(biglmfit)))
#! }
#! \keyword{ IO }
#! \keyword{ data }


chunk.ff_vector <- function(x, RECORDBYTES = .rambytes[vmode(x)], BATCHBYTES = getOption("ffbatchbytes"), ...){
  n <- length(x)
  if (n){
    l <- list(...)
    if (is.null(l$from))
      l$from <- 1L
    if (is.null(l$to))
      l$to <- n
    if (is.null(l$by) && is.null(l$len)){
      b <- BATCHBYTES %/% RECORDBYTES
      if (b==0L){
        b <- 1L
        warning("single record does not fit into BATCHBYTES")
      }
      l$by <- b
    }
    l$maxindex <- n
    ret <- do.call("chunk.default", l)

  }else{
    ret <- list()
  }
  ret
}


chunk.ffdf <- function(x, RECORDBYTES = sum(.rambytes[vmode(x)]), BATCHBYTES = getOption("ffbatchbytes"), ...){
  n <- nrow(x)
  if (n){
    l <- list(...)
    if (is.null(l$from))
      l$from <- 1L
    if (is.null(l$to))
      l$to <- n
    if (is.null(l$by) && is.null(l$len)){
      b <- BATCHBYTES %/% RECORDBYTES
      if (b==0L){
        b <- 1L
        warning("single record does not fit into BATCHBYTES")
      }
      l$by <- b
    }
    l$maxindex <- n
    ret <- do.call("chunk.default", l)

  }else{
    ret <- list()
  }
  ret
}




vw.ffdf <- function(x, ...){
  warning("vw(ffdf) not allowed")
  NULL
}

"vw<-.ffdf" <- function(x, ..., value){
  stop("vw(ffdf)<- not allowed")
}



dim.ffdf <- function(x){
  attr(.subset2(x,"virtual"), "Dim")
}

"dim<-.ffdf" <- function(x, value){
  if (is.null(value) || length(value)!=2)
    stop("wrong value assigned")
  d <- dim(x)
  d1 <- d[[1]]
  d2 <- d[[2]]
  value <- as.integer(value)
  v1 <- value[[1]]
  v2 <- value[[2]]
  if (is.na(v2) && v2!=d2)
    stop("you may only change the number of rows")
  if (v1==d1)
    return(x)

  rnam <- row.names(x)

  physical <- .subset2(x, "physical")


  np <- length(physical)
  for (p in seq_len(np)){
    d <- dim(physical[[p]])
    if (is.null(d))
      length(physical[[p]]) <- v1
    else
      dim(physical[[p]]) <- c(v1, d[[2]])
  }

  oldclass <- oldClass(x)
  oldClass(x) <- NULL
  x$physical <- physical
  attr(x$virtual, "Dim") <- value

  oldClass(x) <- oldclass

  if (!is.null(rnam)){
    n <- length(rnam)
    length(rnam) <- v1
    if (v1>n)
      rnam[(n+1L):v1] <- (n+1L):v1
    row.names(x) <- rnam
  }

  x
}


dimorder.ffdf <- function(x, ...){
  attr(.subset2(x,"virtual"), "Dimorder")
}

# assigning dimorder has no immediate consequences, but might signal to R.ff methods whether prefer to process in columns (1:2) or rows (2:1)
"dimorder<-.ffdf" <- function(x, ..., value){
  value <- as.integer(value)
  if (!identical(sort(value), 1:2))
    stop("illegal dimorder")
  cl <- oldClass(x)
  oldClass(x) <- NULL
  attr(x$virtual, "Dimorder") <- value
  oldClass(x) <- cl
  x
}


#! \name{dimnames.ffdf}
#! \Rdversion{1.1}
#! \alias{dimnames.ffdf}
#! \alias{dimnames<-.ffdf}
#! \alias{names.ffdf}
#! \alias{names<-.ffdf}
#! \alias{row.names.ffdf}
#! \alias{row.names<-.ffdf}
#! \title{
#!   Getting and setting dimnames of ffdf
#! }
#! \description{
#!   Getting and setting dimnames, columnnames or rownames
#! }
#! \usage{
#!   \method{dimnames}{ffdf}(x)
#!   \method{dimnames}{ffdf}(x) <- value
#!   \method{names}{ffdf}(x)
#!   \method{names}{ffdf}(x) <- value
#!   \method{row.names}{ffdf}(x)
#!   \method{row.names}{ffdf}(x) <- value
#! }
#! \arguments{
#!   \item{x}{ a \code{\link{ffdf}} object }
#!   \item{value}{ a character vector, or, for dimnames a list with two character vectors }
#! }
#! \details{
#!   It is recommended not to assign row.names to a large ffdf object.
#! }
#! \value{
#!   The assignment function return the changed ffdf object. The other functions return the expected.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{ffdf}}, \code{\link{dimnames.ff}}, \code{\link{rownames}}, \code{\link{colnames}}
#! }
#! \examples{
#!   ffd <- as.ffdf(data.frame(a=1:26, b=letters))
#!   dimnames(ffd)
#!   row.names(ffd) <- letters
#!   dimnames(ffd)
#!   ffd
#!   rm(ffd); gc()
#! }

names.ffdf <- function(x){
  rownames(.subset2(x,"virtual"))
}

"names<-.ffdf" <- function(x, value){
  cl <- oldClass(x)
  oldClass(x) <- NULL
  rownames(x$virtual) <- value
  oldClass(x) <- cl
  x
}


row.names.ffdf <- function(x){
  .subset2(x,"row.names")
}

"row.names<-.ffdf" <- function(x, value){
  cl <- oldClass(x)
  oldClass(x) <- NULL
  x$row.names <- value
  oldClass(x) <- cl
  x
}

dimnames.ffdf <- function(x){
  list(
    .subset2(x,"row.names")
  , rownames(.subset2(x,"virtual"))
  )
}

"dimnames<-.ffdf" <- function(x, value){
  cl <- oldClass(x)
  oldClass(x) <- NULL
  x$row.names <- value[[1]]
  rownames(x$virtual) <- value[[2]]
  oldClass(x) <- cl
  x
}

#! \name{physical.ffdf}
#! \Rdversion{1.1}
#! \alias{physical.ffdf}
#! \alias{virtual.ffdf}
#! \title{
#! Getting physical and virtual attributes of ffdf objects
#! }
#! \description{
#! Functions for getting physical and virtual attributes of ffdf objects.
#! }
#! \usage{
#! \method{physical}{ffdf}(x)
#! \method{virtual}{ffdf}(x)
#! }
#! \arguments{
#!   \item{x}{an \code{\link{ffdf}} object}
#! }
#! \details{
#! \code{\link{ffdf}} objects enjoy a complete decoupling of virtual behaviour from physical storage.
#! The physical component is simply a (potentially named) list where each element represents an atomic ff vector or matrix.
#! The virtual component is itself a dataframe, each row of which defines a column of the ffdf through a mapping to the physical component.
#! }
#! \value{
#! 'physical.ffdf' returns a \code{\link{list}} with atomic ff objects. \cr
#! 'virtual.ffdf' returns a \code{\link{data.frame}} with the following columns \cr
#!   \item{VirtualVmode}{the \code{\link{vmode}} of this row (=ffdf column)}
#!   \item{AsIs}{logical defining the \code{\link{AsIs}} status of this row (=ffdf column)}
#!   \item{VirtualIsMatrix}{logical defining whether this row (=ffdf column) represents a matrix}
#!   \item{PhysicalIsMatrix}{logical reporting whether the corresponding physical element is a matrix}
#!   \item{PhysicalElementNo}{integer identifying the corresponding physical element}
#!   \item{PhysicalFirstCol}{integer identifying the first column of the corresponding physical element (1 if it is not a matrix)}
#!   \item{PhysicalLastCol}{integer identifying the last column of the corresponding physical element (1 if it is not a matrix)}
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{ffdf}}, \code{\link[=physical.ff]{physical}}, \code{\link[=virtual.ff]{virtual}}, \code{\link[=vmode.ffdf]{vmode}}
#! }
#! \examples{
#!   x <- 1:2
#!   y <- matrix(1:4, 2, 2)
#!   z <- matrix(1:4, 2, 2)
#!
#!   message("Here the y matrix is first converted to single columns by data.frame, 
#! then those columns become ff")
#!   d <- as.ffdf(data.frame(x=x, y=y, z=I(z)))
#!   physical(d)
#!   virtual(d)
#!
#!   message("Here the y matrix is first converted to ff, and then stored still as matrix 
#! in the ffdf object (although virtually treated as columns of ffdf)")
#!   d <- ffdf(x=as.ff(x), y=as.ff(y), z=I(as.ff(z)))
#!   physical(d)
#!   virtual(d)
#!
#!   message("Apply the usual methods extracting physical attributes")
#!   lapply(physical(d), filename)
#!   lapply(physical(d), vmode)
#!   message("And don't confuse with virtual vmode")
#!   vmode(d)
#!
#!   rm(d); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

physical.ffdf <- function(x)
  .subset2(x, "physical")

"physical<-.ffdf" <- function(x, value)
  stop("physical(ffdf) <- not allowed")

virtual.ffdf <- function(x)
  .subset2(x, "virtual")

"virtual<-.ffdf" <- function(x, value)
  stop("virtual(ffdf) <- not allowed")




is.open.ffdf <- function(x, ...){
  p <- .subset2(x, "physical")
  o <- sapply(p, is.open)
  rnam <- .subset2(x, "row.names")
  if (is.ff(rnam))
    o <- c(row.names=is.open(rnam), o)
  if (all(o)) TRUE
  else if (any(o)) NA
  else FALSE
}

open.ffdf <- function (con, readonly = FALSE, pagesize = NULL, caching = NULL
, assert = FALSE
, ...){
  p <- .subset2(con, "physical")
  o <- sapply(p, open, readonly=readonly, pagesize=pagesize, caching=caching, assert=assert)
  rnam <- .subset2(con, "row.names")
  if (is.ff(rnam))
    o <- c(row.names=open(rnam, readonly=readonly, pagesize=pagesize, caching=caching), o)
  if (any(is.na(o))) NA
  else if (all(o)) TRUE
  else if (any(o)) NA
  else FALSE
}

close.ffdf <- function (con, ...){
  p <- .subset2(con, "physical")
  o <- sapply(p, close)
  rnam <- .subset2(con, "row.names")
  if (is.ff(rnam))
    o <- c(row.names=close(rnam), o)
  if (any(is.na(o))) NA
  else if (all(o)) TRUE
  else if (any(o)) NA
  else FALSE
}

delete.ffdf <- function(x, ...){
  p <- .subset2(x, "physical")
  o <- sapply(p, delete)
  rnam <- .subset2(x, "row.names")
  if (is.ff(rnam))
    o <- c(row.names=delete(rnam), o)
  if (any(is.na(o))) NA
  else if (all(o)) TRUE
  else if (any(o)) NA
  else FALSE
}


print.ffdf <- function(x, maxdim=c(16, 16), digits = getOption("digits"), ...){
  v <- .subset2(x, "virtual")
  p <- .subset2(x, "physical")
  rnam <- .subset2(x, "row.names")
  v$PhysicalName <- names(p)[v$PhysicalElementNo]
  v$PhysicalVmode <- sapply(p, vmode)[v$PhysicalElementNo]
  v$PhysicalIsOpen <- sapply(p, is.open)[v$PhysicalElementNo]
  v <- v[,c("PhysicalName", "VirtualVmode", "PhysicalVmode", "AsIs", "VirtualIsMatrix", "PhysicalIsMatrix", "PhysicalElementNo", "PhysicalFirstCol", "PhysicalLastCol", "PhysicalIsOpen")]
  o <- is.open(x)
  cat("ffdf (", if (is.na(o)) "some open" else if (o) "all open" else "all closed", ") dim=c(", paste(dim(x), collapse=","), "), dimorder=c(", paste(dimorder(x), collapse=","), ") row.names=", sep="")
  if (is.null(rnam))
    cat("NULL\n")
  else if (is.ff(rnam)){
    print(rnam)
  }else
    cat(paste(class(rnam), collapse=","), "\n")
  cat("ffdf virtual mapping\n")
  print(v)
  cat("ffdf data\n")
  print(matprint(x, maxdim=maxdim, digits=digits))
  invisible(v)
}

str.ffdf <- function(object, nest.lev=0, ...){
  nest.str <- paste(rep(" ..", nest.lev), collapse="")
  cat(nest.str, "List of 3\n", sep="")
  cat(nest.str, " $ virtual: ", sep="")
  str(.subset2(object, "virtual"), nest.lev=nest.lev+1)
  cat(nest.str, " $ physical: ", sep="")
  str(.subset2(object, "physical"), nest.lev=nest.lev+1)
  cat(nest.str, " $ row.names: ", sep="")
  str(.subset2(object, "row.names"), nest.lev=nest.lev+1)
  cat(nest.str, "- attributes: ", sep="")
  str(attributes(object), nest.lev=nest.lev+1)
}
