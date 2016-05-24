#' Sparsity of Matrix Objects
#' 
#' Determine the sparsity of a matrix, distributed, dense, or otherwise.
#' 
#' The sparsity count of a matrix is returned.
#' 
#' @aliases sparsity-method sparsity,vector-method sparsity,matrix-method
#' sparsity,dmat-method sparsity
#' 
#' @param x 
#' numeric matrix
#' @param count 
#' character; options are "zero" and "other". The former counts
#' the number of zeros, while the latter counts the number of non-zeros
#' ('other' elements).
#' @param out 
#' character; options are "count", "proportion", and "percent". This
#' determines whether a pure count, proportion of \code{count} elements in the
#' matrix, or percentage of \code{count} elements in the matrix.
#' @param tol 
#' numeric; the tolerance for numerical zero. This is ignored if the
#' input data is integer/logical.
#' 
#' @keywords Methods,Sparse
#' 
#' @name sparsity
#' @rdname sparsity
NULL



sparse_count_zeros <- function(x, tol=.Machine$double.eps)
{
  if (is.logical(x))
    storage.mode(x) <- "integer"
  else if (!is.integer(x) && !is.double(x))
    storage.mode(x) <- "double"
  
  if (is.integer(x))
    ret <- .Call(R_int_sparse_count_zeros, x)
  else
    ret <- .Call(R_sparse_count_zeros, x, tol)
  
  return( ret )
}



check_sparsity_inputs <- function(count, out, tol)
{
  pbdMPI::comm.match.arg(tolower(count), c("zero", "other"))
  pbdMPI::comm.match.arg(tolower(out), c("count", "proportion", "percent"))
}



calc_sparsity_return <- function(n, dim, count, out)
{
  count <- pbdMPI::comm.match.arg(tolower(count), c("zero", "other"))
  out <- pbdMPI::comm.match.arg(tolower(out), c("count", "proportion", "percent"))
  
  if (count == "other")
    n <- dim - n
  
  if (out == "count")
    ret <- n
  else if (out == "proportion")
    ret <- n/dim
  else if (out == "percent")
    ret <- n/dim*100
  
  return( ret )
}



#' @rdname sparsity
#' @export
setMethod("sparsity", signature(x="matrix"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    check_sparsity_inputs(count=count, out=out, tol=tol)
    n <- sparse_count_zeros(x=x, tol=tol)
    
    dim <- prod(dim(x))
    ret <- calc_sparsity_return(n=n, dim=dim, count=count, out=out)
    
    return( ret )
  }
)



#' @rdname sparsity
#' @export
setMethod("sparsity", signature(x="vector"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    if (!is.numeric(x))
      comm.stop("argument 'x' must be a numeric vector")
    
    dim(x) <- c(length(x), 1L)
    ret <- sparsity(x=x, count=count, out=out, tol=tol)
    
    return( ret )
  }
)



#' @rdname sparsity
#' @export
setMethod("sparsity", signature(x="dmat"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    n <- sparsity(x=x@Data, count="zero", out="count", tol=tol)
    
    n <- pbdMPI::allreduce(n)
    
    dim <- prod(x@dim)
    ret <- calc_sparsity_return(n=n, dim=dim, count=count, out=out)
    
    return( ret )
  }
)

