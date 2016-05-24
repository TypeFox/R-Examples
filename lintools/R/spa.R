


#' Project a vector on the border of the region defined by a set of linear (in)equality restrictions.
#'
#' Compute a vector, closest to \code{x} in the Euclidean sense, satisfying a set of linear (in)equality restrictions.
#'
#' @param x [\code{numeric}] Vector that needs to satisfy the linear restrictions.
#' @param A [\code{matrix}] Coefficient matrix for linear restrictions.
#' @param b [\code{numeric}] Right hand side of linear restrictions.
#' @param w [\code{numeric}] Optional weight vector of the same length as \code{x}. Must be positive.
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and \code{b} are treated as linear equalities. 
#'    The others as Linear inequalities of the form \eqn{Ax<=b}.
#' @param eps The maximum allowed deviation from the constraints (see details).
#' @param maxiter maximum number of iterations
#'
#' @section Details:
#'
#' The tolerance \code{eps} is defined as the maximum absolute value of the difference vector 
#' \eqn{\boldsymbol{Ax}-\boldsymbol{b}} for equalities. For inequalities, the difference vector
#' is set to zero when it's value is lesser than zero (i.e. when the restriction is satisfied). The
#' algorithm iterates until either the tolerance is met, the number of allowed iterations is
#' exceeded or divergence is detected. 
#' 
#' @return
#' A \code{list} with the following entries:
#' \itemize{
#'  \item{\code{x}: the adjusted vector}
#'  \item{\code{status}: Exit status:
#'   \itemize{
#'    \item{0: success}
#'    \item{1: could not allocate enough memory (space for approximately \eqn{2(m+n)} \code{double}s is necessary).}
#'    \item{2: divergence detected (set of restrictions may be contradictory)}
#'    \item{3: maximum number of iterations reached}
#'   }
#'  }
#'  \item{\code{eps}: The tolerance achieved after optimizing (see Details).}
#'  \item{\code{iterations}: The number of iterations performed.}
#'  \item{\code{duration}: the time it took to compute the adjusted vector}
#'  \item{\code{objective}: The (weighted) Euclidean distance between the initial and the adjusted vector}
#' }
#' @example ../examples/project.R
#' 
#' @seealso \code{\link{sparse_project}}
#' @export
project <- function(x,A,b, neq=length(b), w=rep(1.0,length(x)), eps=1e-2, maxiter=1000L){
  
  check_sys(A=A, b=b, neq=neq, x=x, eps=eps)

  stopifnot(is.numeric(x)
  , length(w) == length(x)
  , is.numeric(w)
  , all_finite(w)
  , all(w>0)
  , maxiter > 0
  , is.finite(maxiter)
  )

  storage.mode(x) <- "double"
   
  t0 <- proc.time()
  y <- .Call("R_dc_solve", 
    A, 
    as.double(b), 
    as.double(w),
    as.integer(neq),
    as.double(eps),
    as.integer(maxiter),
    as.double(x)
  )
  
  t1 <- proc.time()
  objective <- sqrt(sum(w*(x-as.vector(y))^2))
  eps <- attr(y,"eps")
  status <- attr(y,"status")
  niter  <- attr(y,"niter")
  attributes(y) <- NULL
  
  list(x = y
    , status = status
    , eps=eps
    , iterations = niter
    , duration=t1-t0 
    , objective=objective
  )
} 

#' Successive projections with sparsely defined restrictions
#'
#' Compute a vector, closest to \eqn{x} satisfying a set of linear (in)equality restrictions.
#' 
#' 
#' @param x \code{[numeric]} Vector to optimize, starting point.
#' @param A \code{[data.frame]} Coeffiencient matrix in \code{[row,column,coefficient]} format.
#' @param b \code{[numeric]} Constant vector of the system \eqn{Ax\leq b}
#' @param neq \code{[integer]} Number of equalities
#' @param w \code{[numeric]} weight vector of same length of \code{x}
#' @param eps maximally allowed tolerance
#' @param maxiter maximally allowed number of iterations.
#' @param ... extra parameters passed to \code{\link{sparse_constraints}}
#'
#' @section Details:
#'
#' The tolerance \code{eps} is defined as the maximum absolute value of the difference vector 
#' \eqn{\boldsymbol{Ax}-\boldsymbol{b}} for equalities. For inequalities, the difference vector
#' is set to zero when it's value is lesser than zero (i.e. when the restriction is satisfied). The
#' algorithm iterates until either the tolerance is met, the number of allowed iterations is
#' exceeded or divergence is detected. 
#' 
#' @return
#' A \code{list} with the following entries:
#' \itemize{
#'  \item{\code{x}: the adjusted vector}
#'  \item{\code{status}: Exit status:
#'   \itemize{
#'    \item{0: success}
#'    \item{1: could not allocate enough memory (space for approximately \eqn{2(m+n)} \code{double}s is necessary).}
#'    \item{2: divergence detected (set of restrictions may be contradictory)}
#'    \item{3: maximum number of iterations reached}
#'   }
#'  }
#'  \item{\code{eps}: The tolerance achieved after optimizing (see Details).}
#'  \item{\code{iterations}: The number of iterations performed.}
#'  \item{\code{duration}: the time it took to compute the adjusted vector}
#'  \item{\code{objective}: The (weighted) Euclidean distance between the initial and the adjusted vector}
#' }
#' @seealso \code{\link{project}}, \code{\link{sparse_constraints}}
#'
#' @example ../examples/sparse_project.R
#' @export
sparse_project <- function(x, A, b, neq=length(b)
    , w=rep(1.0,length(x)), eps=1e-2, maxiter=1000L,...){
  sc <- sparse_constraints(object=A,b=b,neq=neq,...)
  sc$project(x=x, w=w, eps=eps, maxiter = maxiter)
}


