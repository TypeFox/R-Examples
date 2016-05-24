#' Extract or Replace Parts of a Distributed Matrix
#' 
#' Operators to extract or replace parts of a distributed matrix.
#' 
#' \code{[} can be used to extract/replace for a distributed matrix exactly as
#' you would with an ordinary matrix.
#' 
#' The functions rely on reblocking across different BLACS contexts.  If
#' \code{i} is not empty, then the input distributed matrix will be
#' redistributed along context 1, where extracting/deleting rows does not
#' destroy block-cyclicality. Likewise, if \code{j} is not empty, then the
#' input distributed matrix will be redistributed along context 2. When
#' extraction is complete, the matrix will be redistributed across its input
#' context.
#' 
#' @param x 
#' numeric distributed matrix.
#' @param i,j 
#' indices specifying elements to extract or replace.  Indices can
#' be \code{numeric}, \code{character}, empty, or \code{NULL}.
#' number of elements for a vector (including lists), rows for a matrix or data
#' frame or lines for a function. If negative, all but the \code{n} last/first
#' number of elements of \code{x}.
#' @param ... 
#' additional arguments.
#' @param ICTXT 
#' optional BLACS context number for output
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix(1:9, 3, bldim=2)
#' 
#' y <- x[, -1]
#' y <- head(y, 2)
#' y
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction
#' @name extract
#' @rdname extract
NULL



#' @rdname extract
#' @export
setMethod("[", signature(x="ddmatrix"),
  function(x, i, j, ICTXT)
  {
    attributes(x@Data) <- attributes(x@Data)[which(names(attributes(x@Data))=='dim')]
    
    if (missing(ICTXT))
      oldCTXT <- x@ICTXT
    else
      oldCTXT <- ICTXT
    oldbldim <- x@bldim
    if (missing(i) && missing(j))
      return(x)
    else
      newObj <- x
    
    imiss <- missing(i)
    if (!imiss){
      if (is.ddmatrix(i)){
        if (comm.any(is.logical(i@Data))){
          i <- as.vector(i)
          storage.mode(i) <- "logical"
        }
        else
          i <- as.vector(i)
      }
      if (is.logical(i))
        i <- which(as.vector(i > 0))
      
      ilng <- length(i)
    }
    else
      ilng <- x@dim[1L]
    
    jmiss <- missing(j)
    if (!jmiss){
      if (is.ddmatrix(j)){
        if (comm.any(is.logical(j@Data))){
          j <- as.vector(j)
          storage.mode(j) <- "logical"
        }
        else
          j <- as.vector(j)
      }
      if (is.logical(j))
        j <- which(as.vector(j > 0))
      
      jlng <- length(j)
    }
    else
      jlng <- x@dim[2L]
    
    # special cases 
    if (!imiss && !jmiss){
      # user wants exactly 1 value
      if (ilng==1 && i>0 && jlng==1 && j>0){
        coords <- base.g2l_coord(ind=c(i, j), dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
        if (all(!is.na(coords[c(3,4)])))
          out <- x@Data[coords[5], coords[6]]
        else
          out <- 0
        out <- reduce(out, op='sum')
        if (comm.rank() > 0)
          out <- 0
        out <- new("ddmatrix", Data=matrix(out), dim=c(1,1), ldim=c(1,1), bldim=x@bldim, ICTXT=x@ICTXT)
        return( out )
      }
#      else if (ilng==1){
#        
#      }
#      else if (jlng==1){
#        
#      }
    }
    
    # special cases:  contiguous blocks starting from row/col 1
    if (imiss || ( ilng==length(i) && all(1:ilng == i) ))
    {
      if (jmiss || ( jlng==length(j) && all(1:jlng == j)) )
      {
        # user wants block [1:i] x [1:j]
        dim <- c(ilng, jlng)
        ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT, fixme=TRUE)
        if ( base.ownany(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT) )
          Data <- x@Data[1L:ldim[1L], 1L:ldim[2L], drop=FALSE]
        else 
          Data <- matrix(0, 1, 1)
        
        out <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=x@bldim, ICTXT=x@ICTXT)
        
        return( out )
      }
    }
    
    
    # general case
    if (!imiss) { # skip if no 'i' was supplied
      if (ilng > 0) # ignore i = numeric(0)
#        if (newObj@ICTXT != 1)
          newObj <- base.dropper(x=newObj, oldbldim=oldbldim, iorj='i', ij=i, ICTXT=1)
    }
    
    if (!jmiss){
      if (jlng > 0)
        if (base::length(j)>0)
#          if (newObj@ICTXT != 2)
          newObj <- base.dropper(x=newObj, oldbldim=oldbldim, iorj='j', ij=j, ICTXT=2)
    }
    
    # bring everything back to full process grid
    if (newObj@ICTXT != oldCTXT)
      newObj <- dmat.reblock(dx=newObj, bldim=oldbldim, ICTXT=oldCTXT)
    
    return(newObj)
  }
)



#' Directly Insert Into Distributed Matrix Submatrix Slot
#' 
#' Allows you to directly replace the submatrix of a distributed matrix.
#' 
#' \code{[<-} allows the user to insert values into a distributed matrix in
#' exactly the same way one would with an ordinary matrix. The indices here are
#' global, meaning that \code{x[i, j]} refers to the \code{(i, j)}'th element
#' of the "full", global matrix, and not necessarily the \code{(i, j)}'th
#' element of the local submatrix.
#' 
#' On the other hand, \code{submatrix<-} is different. It is basically
#' syntactic sugar for:
#' 
#' \code{x@@Data <- newMatrix}
#' 
#' It does not alter the distributed matrix \code{x}'s \code{dim} or
#' \code{bldim}. It \emph{does} adjust the \code{ldim} automatically.  However,
#' using this can be dangerous. It is merely provided to give consistent
#' behavior with the \code{submatrix()} function.
#' 
#' @param x 
#' numeric distributed matrix.
#' @param i,j 
#' global integer indices.
#' @param ...
#' Additional arguments.
#' @param value 
#' replacement value. Can be a global vector or a \code{ddmatrix}.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix(1:9, 3, bldim=2)
#' 
#' x[1, ] <- 0
#' comm.print(submatrix(x), all.rank=T)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction
#' @name insert
#' @rdname insert
NULL



#' @rdname insert
#' @export
setReplaceMethod("[", signature(x ="ddmatrix", value="ANY"),
  function(x, i, j, ..., value) 
  {
    if (missing(i))
      i <- 1L:x@dim[1L]
    if (missing(j))
      j <- 1L:x@dim[2L]
    
    if (any(i > x@dim[1L]) || any(j > x@dim[2L]))
      comm.stop("Error : subscript out of bounds")
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    out <- base.rl2insert(x=x@Data, descx=descx, vec=value, i=i, j=j)
    
    x@Data <- out
    
    return(x)
  }
)



#' @rdname insert
#' @export
setReplaceMethod("[", signature(x ="ddmatrix", value="ddmatrix"),
  function(x, i, j, ..., value) 
  {
    if (missing(i) && missing(j))
      return(x)
    if (missing(i)){
      lv <- as.integer(value@dim[2L])
      if (length(j) %% lv != 0)
        comm.stop("number of items to replace is not a multiple of replacement length")
      else if (any(j > x@dim[2L]))
        comm.stop("subscript out of bounds")
      else {
        descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
        descy <- base.descinit(dim=value@dim, bldim=value@bldim, ldim=value@ldim, ICTXT=value@ICTXT)
        
        out <- base.rcolcpy2(x=x@Data, descx=descx, y=value@Data, descy=descy, xcol=j, ycol=1L:lv)
        ret <- x
        ret@Data <- out
      }
    }
    else if (missing(j))
    {
      lv <- as.integer(value@dim[1L])
      if (length(i) %% lv != 0)
        comm.stop("number of items to replace is not a multiple of replacement length")
      else if (any(i > x@dim[1L]))
        comm.stop("subscript out of bounds")
      else {
        descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
        descy <- base.descinit(dim=value@dim, bldim=value@bldim, ldim=value@ldim, ICTXT=value@ICTXT)
        
        out <- base.rrowcpy2(x=x@Data, descx=descx, y=value@Data, descy=descy, xrow=i, yrow=1L:lv)
        ret <- x
        ret@Data <- out
      }
    }
    
    return( ret )
  }
)

