#' Accessor Functions for Distributed Matrix Slots
#' 
#' Functions to get dimension information, local storage, or current BLACS
#' context from a distributed matrix.
#' 
#' The functions \code{nrow()}, \code{ncol()}, \code{length()} and \code{dim()}
#' are the natural extensions of their ordinary matrix counterparts.
#' 
#' \code{ldim()} will give the dimension of the matrix stored locally on the
#' process which runs the function. This is a local value, so its return is
#' process-dependent.  For example, if the 3x3 global matrix \code{x} is
#' distributed as the \code{ddmatrix} \code{dx} across two processors with
#' process 0 owning the first two rows and process 1 owning the third, then
#' \code{ldim(dx)} will return \code{2 3} on process 0 and \code{1 3} on
#' process 1.
#' 
#' \code{bldim()} will give the blocking dimension that was used to
#' block-cyclically distribute the distributed matrix.
#' 
#' \code{submatrix()} will give the local storage for the requested object.
#' 
#' \code{ICTXT()} will give the current BLACS context (slot ICTXT) for the
#' requested object.
#' 
#' \code{ownany()} is intended mostly for developers.  It answers the question
#' "do I own any of the data?".  The user can either pass a distributed matrix
#' object or the dim, bldim, and ICTXT of one.
#' 
#' @param x 
#' numeric distributed matrix
#' @param dim 
#' global dimension.
#' @param bldim 
#' blocking dimension.
#' @param ICTXT 
#' BLACS context.
#' @param ... 
#' Extra arguments.
#' 
#' @return 
#' Each of \code{dim()}, \code{ldim()}, \code{bldim()} return a length
#' 2 vector.
#' 
#' Each of \code{nrow()}, \code{ncol()}, and \code{length()} return a length 1
#' vector. Likewise, so does \code{ICTXT()}.
#' 
#' \code{submatrix()} returns a matrix; namely, \code{submatrix(x)} returns a
#' matrix of dimensions \code{ldim(x)}.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
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
#' y <- list(dim=dim(x), ldim=ldim(x), bldim=bldim(x))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name Accessors
#' @rdname accessors
NULL

#' @rdname accessors
#' @export
setGeneric(name = "nrow", useAsDefault = base::nrow, package="pbdDMAT")

#' @rdname accessors
#' @export
setMethod("nrow", signature(x="ddmatrix"),
  function(x)
    return(x@dim[1L])
)

#' @rdname accessors
#' @export
setGeneric(name = "NROW", useAsDefault = base::NROW, package="pbdDMAT")

#' @rdname accessors
#' @export
setMethod("NROW", signature(x="ddmatrix"),
  function(x)
    return(x@dim[1L])
)

#' @rdname accessors
#' @export
setGeneric(name = "ncol", useAsDefault = base::ncol, package="pbdDMAT")

#' @rdname accessors
#' @export
setMethod("ncol", signature(x="ddmatrix"),
  function(x)
    return(x@dim[2L])
)

#' @rdname accessors
#' @export
setGeneric(name = "NCOL", useAsDefault = base::NCOL, package="pbdDMAT")

#' @rdname accessors
#' @export
setMethod("NCOL", signature(x="ddmatrix"),
  function(x)
    return(x@dim[2L])
)

#' @rdname accessors
#' @export
setGeneric(name="submatrix", 
  function(x, ...) 
    standardGeneric("submatrix"), 
  package="pbdDMAT"
)

#' @rdname accessors
#' @export
setMethod("submatrix", signature(x="ddmatrix"),
  function(x)
  {
    if (!is.ddmatrix(x))
      comm.stop("Not a distributed matrix")
    else
      return(x@Data)
  }
)

#' @rdname accessors
#' @export
setGeneric(name="ldim", 
  function(x, ...) 
    standardGeneric("ldim"), 
  package="pbdDMAT"
)

#' @rdname accessors
#' @export
setMethod("ldim", signature(x="ddmatrix"),
  function(x)
  {
    if (!is.ddmatrix(x))
      comm.stop("Not a distributed matrix")
    else
      return(x@ldim)
  }
)

#' @rdname accessors
#' @export
setGeneric(name="bldim", 
  function(x, ...) 
    standardGeneric("bldim"), 
  package="pbdDMAT"
)

#' @rdname accessors
#' @export
setMethod("bldim", signature(x="ddmatrix"),
  function(x)
  {
    if (!is.ddmatrix(x))# && !is.ddvector(x))
      comm.stop("'bldim' only applies objects of class 'ddmatrix' and 'ddvector'")
    else
      return(x@bldim)
  }
)

#' @rdname accessors
#' @export
setGeneric(name="ICTXT", 
  function(x, ...) 
    standardGeneric("ICTXT"), 
  package="pbdDMAT"
)

#' @rdname accessors
#' @export
setMethod("ICTXT", signature(x="ddmatrix"),
  function(x)
  {
    if (!is.ddmatrix(x))# && !is.ddvector(x))
      comm.stop("'ICTXT' only applies objects of class 'ddmatrix' and 'ddvector'")
    else
      return(x@ICTXT)
  }
)



#' @rdname accessors
#' @export
setMethod("dim", signature(x="ddmatrix"),
  function(x)
    return(x@dim)
)

#' @rdname accessors
#' @export
setMethod("length", signature(x="ddmatrix"),
  function(x)
    return(prod(x@dim))
)


