#' Arithmetic Reductions: Sums, Means, and Prods
#' 
#' Arithmetic reductions for distributed matrices.
#' 
#' Performs the reduction operation on a distributed matrix.
#' 
#' There are several legitimately new operations, including \code{rowMin()},
#' \code{rowMax()}, \code{colMin()}, and \code{colMax()}.  These
#' implementations are not really necessary in R because one can easily (and
#' reasonably efficiently) do something like
#' 
#' \code{apply(X=x, MARGIN=1L, FUN=min, na.rm=TRUE)}
#' 
#' But \code{apply()} on a \code{ddmatrix} is \emph{very} costly, and should be
#' used sparingly.
#' 
#' \code{sd()} will compute the standard deviations of the columns, equivalent
#' to calling \code{apply(x, MARGIN=2, FUN=sd)} (which will work for
#' distributed matrices, by the way). However, this should be much faster and
#' use less memory than \code{apply()}.  If \code{reduce=FALSE} then the return
#' is a distributed matrix consisting of one (global) row; otherwise, an
#' \code{R} vector is returned, with ownership of this vector determined by
#' \code{proc.dest}.
#' 
#' @param x 
#' numeric distributed matrix
#' @param na.rm 
#' logical. Should missing (including \code{NaN}) be removed?
#' @param ... 
#' additional arguments
#' 
#' @return 
#' Returns a global numeric vector.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- sum(colMeans(x))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name reductions
#' @rdname reductions
NULL



#' @rdname reductions
#' @export
setGeneric(name="rowMin", 
  function(x, ...)
    standardGeneric("rowMin"),
  package="pbdDMAT"
)

#' @rdname reductions
#' @export
setGeneric(name="rowMax", 
  function(x, ...)
    standardGeneric("rowMax"),
  package="pbdDMAT"
)

#' @rdname reductions
#' @export
setGeneric(name="colMin", 
  function(x, ...)
    standardGeneric("colMin"),
  package="pbdDMAT"
)

#' @rdname reductions
#' @export
setGeneric(name="colMax", 
  function(x, ...)
    standardGeneric("colMax"),
  package="pbdDMAT"
)


# -------------------
# MPI-like BLACS reductions
# -------------------


# Higher level reduction interface
dmat.blacsreduction <- function(x, SCOPE, op, ICTXT, proc.dest=-1, check=TRUE)
{
  if (!is.character(SCOPE)){
    if (SCOPE==1)
      SCOPE <- 'Row'
    else if (SCOPE==2)
      SCOPE <- 'Col'
    else 
      comm.stop("ERROR : invalid argument 'scope'")
  }
  else if (SCOPE != 'Row' && SCOPE != 'Col')
    comm.stop("ERROR : invalid argument 'scope'")
  
  if (is.matrix(x)){
    m <- dim(x)[1L]
    n <- dim(x)[2L]
  }
  else if (is.vector(x)){
    m <- length(x)[1L]
    n <- 1L
  }
  else 
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (length(proc.dest==1)){
    if (proc.dest==-1){
      rdest <- cdest <- -1
    } else {
      proc.dest <- base.pcoord(ICTXT=ICTXT, PNUM=proc.dest)
      rdest <- proc.dest[[1L]]
      cdest <- proc.dest[[2L]]
    }
  }
  else {
    rdest <- proc.dest[1L]
    cdest <- proc.dest[2L]
  }
  
  
  # checking that all m and n agree
  if (check){
    if (SCOPE=='All')
      mx <- pbdMPI::allreduce(c(m,n), op='max')
    else
      mx <- base.igamx2d(ICTXT=ICTXT, SCOPE=SCOPE, m=2L, n=1L, x=c(m,n), lda=2, RDEST=-1, CDEST=-1)
    
    dm <- mx[1L] - m
    dn <- mx[2L] - n
    
    if (dm > 0 || dn > 0){
      dim(x) <- NULL
      
      if (is.integer(x))
        nd <- 0L
      else
        nd <- 0.0
      
      
      x <- c(x, rep(nd, prod(mx)-(m*n)))
      m <- mx[1L]
      n <- mx[2L]
    }
  }
  
  
  if (op == 'sum'){
    if (is.integer(x))
      out <- base.igsum2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
    else
      out <- base.dgsum2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
  }
  else if (op == 'max'){
    if (is.integer(x))
      out <- base.igamx2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
    else
      out <- base.dgamx2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
  }
  else if (op == 'min'){
    if (is.integer(x))
      out <- base.igamn2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
    else
      out <- base.dgamn2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
  }
  else 
    comm.stop("ERROR : invalid argument 'op'")
  
  return( out )
}



dmat.allcolreduce <- function(x, op='sum', ICTXT=.pbd_env$ICTXT)
{
  dmat.blacsreduction(x=x, SCOPE='Col', op=op, ICTXT=ICTXT, proc.dest=-1)
}

allcolreduce <- dmat.allcolreduce


dmat.allrowreduce <- function(x, op='sum', ICTXT=.pbd_env$ICTXT)
{
  dmat.blacsreduction(x=x, SCOPE='Row', op=op, ICTXT=ICTXT, proc.dest=-1)
}

allrowreduce <- dmat.allrowreduce


dmat.colreduce <- function(x, op='sum', proc.dest=0, ICTXT=.pbd_env$ICTXT)
{
  dmat.blacsreduction(x=x, SCOPE='Col', op=op, ICTXT=ICTXT, proc.dest=proc.dest)
}

colreduce <- dmat.colreduce


dmat.rowreduce <- function(x, op='sum', proc.dest=0, ICTXT=.pbd_env$ICTXT)
{
  dmat.blacsreduction(x=x, SCOPE='Row', op=op, ICTXT=ICTXT, proc.dest=proc.dest)
}

rowreduce <- dmat.rowreduce



dmat.rcsum <- function(x, na.rm=FALSE, SCOPE, MEAN=FALSE)
{
  if (SCOPE == 'Row'){
    if (MEAN)
      Data <- rowSums(x@Data / as.double(x@dim[2L]), na.rm=na.rm)
    else
      Data <- rowSums(x@Data, na.rm=na.rm)
    
    dim(Data) <- c(base::length(Data), 1L)
    
    if (x@dim[2L]==1)
      return( Data )
    else
      n <- nrow(Data)
  }
  else {
    if (MEAN)
      Data <- colSums(x@Data / as.double(x@dim[1L]), na.rm=na.rm)
    else
      Data <- colSums(x@Data, na.rm=na.rm)
    
    dim(Data) <- c(1L, base::length(Data))
    
    if (x@dim[1L]==1)
      return( Data )
    else
      n <- ncol(Data)
  }
  
  
  out <- dmat.blacsreduction(x=Data, SCOPE=SCOPE, op='sum', ICTXT=x@ICTXT, proc.dest=-1, check=TRUE)
  
  return( out )
}


dmat.rcminmax <- function(x, na.rm=FALSE, SCOPE, op)
{
  op <- pbdMPI::comm.match.arg(op, c("min", "max"))
  Rop <- eval(parse(text=op))
  
  if (SCOPE == 'Row'){
    Data <- apply(x@Data, 1L, Rop, na.rm=na.rm)
    
    dim(Data) <- c(base::length(Data), 1L)
    
    if (x@dim[2L]==1)
      return( Data )
    else
      n <- nrow(Data)
  }
  else {
    Data <- apply(x@Data, 2L, Rop, na.rm=na.rm)
    
    dim(Data) <- c(1L, base::length(Data))
    
    if (x@dim[1L]==1)
      return( Data )
    else
      n <- ncol(Data)
  }
  
  # have to account for possible lack of ownership
  if (SCOPE=='All')
    mx <- pbdMPI::allreduce(dim(Data), op='max')
  else
    mx <- base.igamx2d(ICTXT=x@ICTXT, SCOPE=SCOPE, m=2L, n=1L, x=dim(Data), lda=2, RDEST=-1, CDEST=-1)
  
  m <- dim(x@Data)[1L]
  n <- dim(x@Data)[2L]
  
  dm <- mx[1L] - m
  dn <- mx[2L] - n
  
  if (!is.double(Data))
    storage.mode(Data) <- "double"
  
  iown <- base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
  
  if (op=='min')
    nd <- pbdMPI::allreduce(max(Data), op='max')
  else if (op=='max')
    nd <- pbdMPI::allreduce(min(Data), op='min')
  
  if (dm > 0 || dn > 0 || !iown){
    dim(Data) <- NULL
    
#    if (op=='min')
#      nd <- Inf
#    else if (op=='max')
#      nd <- -Inf
    
    if(!iown)
      Data <- rep(nd, prod(mx))
    else
      Data <- c(Data, rep(nd, prod(mx)-(m*n)))
    
    dim(Data) <- mx
  }
  
  out <- dmat.blacsreduction(x=Data, SCOPE=SCOPE, op=op, ICTXT=x@ICTXT, proc.dest=-1)
  
  return( out )
}

# -------------------
# Reductions
# -------------------

#' @rdname reductions
#' @export
setMethod("rowSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Row', MEAN=FALSE)
#    dim(Data) <- c(base::length(Data), 1L)
    
    z <- new("ddmatrix", Data=Data, dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("colSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Col', MEAN=FALSE)
    
    z <- new("ddmatrix", Data=Data, dim=c(1L, x@dim[2L]), ldim=c(1L,length(x@Data)), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("rowMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Row', MEAN=TRUE)
#    dim(Data) <- c(base::length(Data), 1L)
    
    z <- new("ddmatrix", Data=Data, dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("colMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Col', MEAN=TRUE)
    
    z <- new("ddmatrix", Data=Data, dim=c(1, x@dim[2]), ldim=c(1,length(x@Data)), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("rowMin", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcminmax(x=x, na.rm=na.rm, SCOPE='Row', op='min')
    
    z <- new("ddmatrix", Data=Data, dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("rowMin", signature(x="matrix"), 
  function(x, na.rm=FALSE)
    apply(X=x, MARGIN=1L, FUN=min, na.rm=na.rm)
)

#' @rdname reductions
#' @export
setMethod("colMin", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcminmax(x=x, na.rm=na.rm, SCOPE='Col', op='min')
    
    z <- new("ddmatrix", Data=Data, dim=c(1L, x@dim[2L]), ldim=c(1L,length(x@Data)), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("colMin", signature(x="matrix"), 
  function(x, na.rm=FALSE)
    apply(X=x, MARGIN=2L, FUN=min, na.rm=na.rm)
)

#' @rdname reductions
#' @export
setMethod("rowMax", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcminmax(x=x, na.rm=na.rm, SCOPE='Row', op='max')
    
    z <- new("ddmatrix", Data=Data, dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("rowMax", signature(x="matrix"), 
  function(x, na.rm=FALSE)
    apply(X=x, MARGIN=1L, FUN=max, na.rm=na.rm)
)

#' @rdname reductions
#' @export
setMethod("colMax", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcminmax(x=x, na.rm=na.rm, SCOPE='Col', op='max')
    
    z <- new("ddmatrix", Data=Data, dim=c(1L, x@dim[2L]), ldim=c(1L,length(x@Data)), bldim=x@bldim) 
    
    return( z )
  }
)

#' @rdname reductions
#' @export
setMethod("colMin", signature(x="matrix"), 
  function(x, na.rm=FALSE)
    apply(X=x, MARGIN=2L, FUN=max, na.rm=na.rm)
)
