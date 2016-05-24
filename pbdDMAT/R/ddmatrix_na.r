#' Handle Missing Values in Distributed Matrices
#' 
#' Dealing with NA's and NaN's.
#' 
#' Removes rows containing NA's and NaN's.
#' 
#' The function relies on reblocking across different BLACS contexts.  The
#' input distributed matrix will be redistributed along context 1, where
#' extracting/deleting rows does not destroy block-cyclicality.
#' 
#' Only advanced users should supply an \code{ICTXT} value. Most should simply
#' leave this argument blank.
#' 
#' The context of the return is dependent on the function arguments.  If the
#' \code{ICTXT=} argument is missing, then the return will be redistributed
#' across its input context \code{object@@ICTXT}.  Otherwise, the return will be
#' redistributed across the supplied \code{ICTXT}.
#' 
#' @param object 
#' numeric distributed matrix
#' @param ... 
#' extra arguments
#' @param ICTXT 
#' optional BLACS context number for output
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
#' x[1, 1] <- NA
#' x <- as.ddmatrix(x)
#' 
#' y <- na.exclude(x)
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction Type
#' @name na
#' @rdname na
NULL



#' @rdname na
#' @export
setGeneric(name = "na.exclude", useAsDefault = stats::na.exclude, package="pbdDMAT")



#' @rdname na
#' @export
setMethod("na.exclude", signature(object="ddmatrix"),
  function(object, ..., ICTXT)
  {
    # 1xn's have to be handled separately
    if (object@dim[1] == 1){
      anynas <- any(is.na(object@Data))
      anynas <- as.logical(allreduce(anynas, op='max'))
      if (anynas){
        object@Data <- matrix(0)
        object@dim[1] <- 0
        object@ldim <- c(1,1)
        if (!missing(ICTXT))
          object@ICTXT <- ICTXT
      } else if (object@ICTXT != ICTXT)
        object <- dmat.reblock(dx=object, bldim=object@bldim, ICTXT=ICTXT)
      
      return(object)
    }
    
    # General case
    if (missing(ICTXT))
      oldCTXT <- object@ICTXT
    else
      oldCTXT <- ICTXT
    blacs_ <- base.blacs(1)

    oldbldim <- object@bldim
    bldim <- c(dim(object)[1], ceiling(oldbldim[2] / blacs_$NPCOL))

    if (object@ICTXT != 1)
      newObj <- dmat.reblock(dx=object, bldim=bldim, ICTXT=1)

    iown <- ownany(dim=newObj@dim, bldim=newObj@bldim, ICTXT=1)

#    if (blacs_$MYROW != -1 && blacs_$MYCOL != -1)   FIXME

    if (iown)
      tmp <- base::rowSums(newObj@Data)
    else
      tmp <- numeric(0)

    tmplen <- pbdMPI::allreduce(length(tmp), op='max')
    if (length(tmp) < tmplen)
      tmp <- rep(0, tmplen)
    tmp <- pbdMPI::allreduce(tmp)

    narows <- which(is.na(tmp))
    lnarows <- length(narows)
    if (lnarows > 0 && iown){
      if (lnarows < newObj@dim[1])
        new <- newObj@Data[-narows, , drop=FALSE] 
      else
        new <- matrix(0.0, nrow=0, ncol=newObj@dim[2])
#        if (!is.matrix(new))
#          new <- matrix(new, nrow=1)
      newObj@Data <- new
      attr(narows, "class") <- "exclude"
      attr(newObj@Data, "na.action") <- narows
    }

    newObj@ldim <- dim(newObj@Data)

    # correction for 0xn ldims
    if (newObj@ldim[1]==0){
      newObj@Data <- matrix(0.0)
      newObj@dim[1] <- 0
      newObj@ldim <- c(1,1)
      newObj@ICTXT <- oldCTXT
    }

    if (all(newObj@dim>0)){
      newdim <- pbdMPI::allreduce(dim(newObj@Data)[1], op='max')
      newObj@dim[1]  <- newdim
    }

    if (newObj@ICTXT != oldCTXT)
      newObj <- dmat.reblock(dx=newObj, bldim=oldbldim, ICTXT=oldCTXT)

    return(newObj)
  }
)

