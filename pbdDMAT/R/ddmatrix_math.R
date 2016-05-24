#' Miscellaneous Mathematical Functions
#' 
#' Binary operations for distributed matrix/distributed matrix and distributed
#' matrix/vector operations.
#' 
#' Performs the miscellaneous mathematical calculation on a distributed matrix.
#' 
#' @param x 
#' numeric distributed matrix
#' @param base 
#' a positive number: the base with respect to which logarithms are
#' computed. Defaults to e='exp(1)'.
#' 
#' @return Returns a distributed matrix.
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
#' y <- sqrt(abs(log(x/10)))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name math
#' @rdname math
NULL



#' @rdname math
#' @export
setMethod("sqrt", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- sqrt(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("abs", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- abs(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("exp", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- exp(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("log", signature(x="ddmatrix"),
  function(x, base=exp(1))
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log(x@Data, base)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("log2", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log2(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("log10", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log10(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("log1p", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log1p(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("sin", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- sin(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("cos", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- cos(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("tan", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- cos(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("asin", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- asin(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("acos", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- acos(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("atan", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- atan(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("sinh", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- sinh(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("cosh", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- cosh(x@Data)
    return(x)
  }
)

#' @rdname math
#' @export
setMethod("tanh", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- tanh(x@Data)
    return(x)
  }
)
