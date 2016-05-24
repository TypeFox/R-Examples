#' @title Nonlinear optimizers using trust regions.
#'
#' @description 
#' Run nonlinear minimizer using trust region algorithm with conjugate
#' gradient search directions and quasi-Hessian updates.
#'
#' @param x A numeric vector of starting values for the optimizer.
#' @param fn An R function that takes \code{x} as its first argument.
#'   Returns the value of the objective function at \code{x}.
#'   Note that the optimizer will \emph{minimize} \code{fn} (see
#'   function.scale.factor under control)
#' @param gr An R function that takes x as its first argument.  Returns a
#'   numeric vector that is the gradient of \code{fn} at \code{x}. The
#'   length of the gradient must be the same as the length of \code{x}.
#'   The user must supply this function.  If an analytic gradient is not
#'   available, and the method is \code{SR1} or \code{BFGS}, the user should consider a
#'   numerical approximation using finite differencing (see the
#'   numDeriv package).  Do not use a finite-differenced gradient with 
#'   the \code{Sparse} method.  That will cause a world of hurt.
#' @param hs An R function that takes x as its first argument.
#'   Returns a Hessian matrix object of class \code{dgCMatrix} (see the \pkg{Matrix} package).
#'   This function is called only if the selected method is \code{Sparse}.    
#' @param method Valid arguments are \code{SR1},\code{BFGS},and \code{Sparse}.
#' @param control A list containing control parameters for the optimizer.
#'   See details.
#' @param ... Additional arguments passed to \code{fn}, \code{gr} and \code{hs}.
#'   All arguments must be named.
#'
#' @return List containing the following items:
#'   \item{fval}{Value of the objective function}
#'   \item{solution}{Parameter vector at the optimum}
#'   \item{gradient}{Gradient at the optimum} 
#'   \item{hessian}{Estimate of the Hessian at the optimum (as class
#'     \code{symmetricMatrix}, returned only for \code{Sparse} method).}
#'   \item{iterations}{Number of iterations before stopping}
#'   \item{status}{A message describing the last state of the iterator}
#'   \item{nnz}{For the Sparse method only, the number of nonzero elements in the lower triangle of the Hessian}.
#' 
#' @section Details:
#'  The following sections explain how to use the package as a whole.
#' 
#' @section Control parameters:
#' The control list should include the following parameters.
#' \describe{
#' \item{start.trust.radius}{Initial radius of the trust region. Default is 5.  If the algorithm returns non-finite values of the objective function early in the process, try a lower number.}
#' \item{stop.trust.radius}{Minimum radius of trust region.  Algorithm will terminate if radius is below this value.  This is because it may not be possible to get the norm of the gradient smaller than prec, and this is another way to get the algorithm to stop.}
#' \item{cg.tol}{tolerance for the conjugate gradient algorithm that is used for the trust region subproblem.  Set it to something very small.  Default is sqrt(.Machine$double.eps)}
#' \item{prec}{Precision for how close the norm of the gradient at the
#'     solution should be to zero, before the algorithm halts.  It is possible that the algorithm
#'     will not get that far, so it will also stop when the radius of the
#'     trust region is smaller thanstop.trust.radius.  If the trust
#'     region radius collapses, but the norm of the gradient really isn't
#'     close to zero, then something terrible has happened.}
#' \item{report.freq}{An integer. The frequency at which the algorithm
#'   will display the current iteration number or function value, among
#'   other things (see \code{report.level}).  Defaults to 1.}
#' \item{report.level}{The amount of detail in each report.  Defaults to 2.}
#' \item{report.precision}{The number of significant digits used in each
#'     report. Defaults to 5.}
#' \item{maxit}{Maximum number of iterations.  Defaults to 100.}
#' \item{contract.factor}{When the algorithm decides to shrink the trust
#'     region, it will multiply the trust radius by this factor. Defaults to 0.5.}
#' \item{expand.factor}{When the algorithm decides to expand the trust
#'     region, it will multiply the algorithm by this factor. Defaults to 3.}
#' \item{contract.threshold}{The algorithm with accept a proposed move if the ratio of the actual improvement in the objective function, to the predicted improvement from the trust region subproblem, is greater than this amount.  Otherwise, the trust region will contract.  Default is 0.25.}
#' \item{expand.threshold.ap}{First criterion to determine if the trust region should expand.  If the ratio of the actual and proposed improvements in the objective function is less than this factor, the algorithm will consider expanding the trust region.  See \code{expand.threshold.radius}. Default is 0.8.}
#' \item{expand.threshold.radius}{If the ratio of the actual and proposed improvement in the objective function is less than \code{expand.threshold.ap}, then, if the normed distance of the proposed move is greater than  \code{expand.threshold.radius}, times the current trust region radius, the trust region will expand.  Default is 0.8.}
#' \item{function.scale.factor}{The algorithm will minimize \code{fn} times this
#'     factor. If you want to maximize \code{fn}, this value should be negative
#'     (usually -1).  Default is 1.}
#' \item{precond.refresh.freq}{Frequency at which the preconditioner
#'     for the conjugate gradiate estimation of the trust region
#'     subproblem is reestimated.  Preconditioners can help the convergence properties of the algorithm.  Default is 1.}
#' \item{preconditioner}{ID for choice of preconditioner.  0 is the
#'     identity matrix (default), For the \code{Sparse} method, 1 is a modified Cholesky preconditioner. For the \code{BFGS} method, 1 is the full Cholesky decomposition.  If you select 1 for the \code{SR1} method, the algorithm will use the identity preconditioner instead.}
#' \item{trust.iter}{Maximum number of conjugate gradient iterations to run when solving the trust region subproblem.  A higher number will lead to more accurate solutions to the subproblem, but may also lead to longer run times. Defaults to 2000.} 
#' }
#'
#' @section Report levels:
#' The \code{report.level} control parameter determines how much information is displayed each time the algorithm reports the current state.  Possible values are 
#' 
#' \describe{
#' \item{<=0}{No information (a quiet run)}
#' \item{1}{Current iteration number, and current value of the objective function.}
#' \item{2}{Information from level 1, plus the current norm of the gradient and a status message.}
#' \item{3}{Information from levels 1 and 2, plus the current normed radius of the trust region.}
#' \item{4}{Information from levels 1, 2, and 3, plus information from each estimate of the trust region subproblem (number of conjugate gradient iterations and how/why the CG algorithm terminated).}
#' }
#' 
#' Default level is 2.  Levels 3 and 4 are available primarily for debugging purposes.
#' 
#' @section Stopping criteria:
#' The algorithm will stop when one of the following three conditions are met:
#' \itemize{
#' \item{The norm of the gradient, divided by the square root of the number of parameters, is less than \code{prec}.}
#' \item{The trust region collapse to a radius smaller than machine precision}
#' \item{The algorithm proposes zero or negative improvement in the objective function (should never happen)}
#' \item{The number of iterations reaches the control parameter \code{maxit}}
#' }
#' 
#' If the algorithm appears to have stopped prematurely (i.e., the norm of the gradient is still too large), then one might just restart the algorithm.  For the quasi-Newton algorithms (\code{SR1} and \code{BFGS}), this will refresh the Hessian, and might allow more progress to be made.
#'
#' @section Estimating a sparse Hessian:
#' Sometimes estimating the Hessian is easy (e.g., you have an analytic representation, or you are using some kind of algorithmic differentiation software).  If you do not know the Hessian, but you do know the sparsity structure, try the \pkg{sparseHessianFD} package. The routines in \pkg{sparseHessianFD} compute the Hessian using finite differencing, but in a way that exploits the sparsity structure.  In many cases, this can be faster than constructing an analytic Hessian for a large problem (e.g., when the Hessian has a block-arrow structure with a large number of blocks).
#' 
#' To use the \pkg{sparseHessianFD} package, you need to provide the row and column indices of the non-zero elements of the lower triangle of the Hessian. This structure cannot change during the course of the trust.optim routine.  Also, you really should provide an analytic gradient.  \pkg{sparseHessianFD} computes finite differences of the gradient, so if the gradient itself is finite-differenced, so much error is propagated through that the Hessians are nearly worthless close to the optimum.
#' 
#' Of course, \pkg{sparseHessianFD} is useful only for the \code{Sparse} method.  That said, one may still get decent performance using these routines even if the Hessian is sparse, if the problem is not too large.  Just treat the Hessian as if it were sparse.
#' @export
trust.optim <- function(x, fn, gr, hs=NULL, method=c("SR1","BFGS","Sparse"), control = list(), ...)
{

  if (is.null(method) || (match(method,c("SR1","BFGS","Sparse"),nomatch=0)==0)) {
    stop("Error in trust.optim:  method must be SR1, BFGS, or Sparse.")
  }
  
  if (!is.function(fn)) stop ("Error in trust.optim:  fn must be a function")
  if (!is.function(gr)) stop ("Error in trust.optim:  gr must be a function")
  
  
  fn1 <- function(x) fn(x,...)  ## currying the data in
  gr1 <- if (!is.null(gr)) function(x) gr(x,...)

    

## test fn and gr

  r1 <- fn1(x)
  if (!is.finite(r1)) stop("Error in trust.optim:  fn at starting values is not finite.")
  r1 <- gr1(x)
  if (any(!is.finite(r1))) stop("Error in trust.optim:  at least one element of gr is not finite.")
  if (length(r1)!=length(x)) stop("Error in trust.optim:  length of gradient does not match length of starting values.")

  if(method=="Sparse") {
    hs1 <- if (!is.null(hs)) function(x) hs(x,...)
    r1 <- hs1(x)
    if (class(r1)!="dgCMatrix") stop("Error in trust.optim:  hs function must return object of class dgCMatrix")
  }
  
    ## Defaults :

  con <- list(start.trust.radius=5.0,
              stop.trust.radius=sqrt(.Machine$double.eps),
              cg.tol=sqrt(.Machine$double.eps),
              prec=sqrt(.Machine$double.eps),
              report.freq=1L,
              report.level=2L,
              report.precision=6L,
              maxit=100L,
              contract.factor = 0.5,
              expand.factor = 3,
              contract.threshold = 0.25,
              expand.threshold.ap = 0.8,
              expand.threshold.radius = 0.8,
              function.scale.factor = as.numeric(1),   
              precond.refresh.freq=1L,
              preconditioner=0L,
              quasi.newton.method=0L,
              trust.iter=2000L
              )

    nmsC <- names(con)

    con[(namc <- names(control))] <- control


  con$report.freq <- as.integer(con$report.freq)
  con$report.level <- as.integer(con$report.level)
  con$report.precision <- as.integer(con$report.precision)
  con$precond.refresh.freq <- as.integer(con$precond.refresh.freq)
  con$maxit <- as.integer(con$maxit)
  
## check control parameter values here     

  if (!is.numeric(con$start.trust.radius) || con$start.trust.radius<=0 || !is.finite(con$start.trust.radius)) {
    stop("Error in trust.optim:  bad value for start.trust.radius.")
  }

  if (!is.numeric(con$stop.trust.radius) || (con$stop.trust.radius <= 0) || !is.finite(con$stop.trust.radius)) {
    stop("Error in trust.optim:  bad value for stop.trust.radius.")
  }

  if (!is.numeric(con$cg.tol) || con$cg.tol<=0 || !is.finite(con$cg.tol)) {
    stop("Error in trust.optim:  bad value for cg.tol.")
  }

  if (!is.numeric(con$prec) || con$prec<=0 || !is.finite(con$prec)) {
    stop("Error in trust.optim:  bad value for prec.")
  }

  if (!is.integer(con$report.freq)) {
    stop("Error in trust.optim:  report.freq must be an integer.")
  }

  if (!is.integer(con$report.level)) {
    stop("Error in trust.optim:  report.level must be an integer.")
  }

  if (!is.integer(con$report.precision)) {
    stop("Error in trust.optim:  report.precision must be an integer.")
  }

  if (!is.integer(con$maxit) || con$maxit<0) {
    stop("Error in trust.optim:  maxit must be an non-negative integer.")
  }

  if (!is.numeric(con$contract.threshold) || con$contract.threshold<0) {
    stop("Error in trust.optim:  contract.threshold must be an non-negative numeric.")
  }

  if (!is.numeric(con$expand.threshold.ap) || con$expand.threshold.ap<0) {
    stop("Error in trust.optim:  expand.threshold.ap must be an non-negative numeric.")
  }

  if (!is.numeric(con$expand.threshold.radius) || con$expand.threshold.radius<0) {
    stop("Error in trust.optim:  expand.threshold.radius must be an non-negative numeric.")
  }

  if (!is.numeric(con$expand.factor) || con$expand.factor<0) {
    stop("Error in trust.optim:  expand.factor must be an non-negative numeric.")
  }

 

  if (method=="Sparse") {
     
    con$quasi.newton.method <- 0L    
##    res <- .Call("sparseTR", x, fn1, gr1, hs1, con)
    res <- sparseTR(x, fn1, gr1, hs1, con)
    res$hessian <- Matrix::t(as(res$hessian,"symmetricMatrix"))
    
  }

  if (method=="SR1" || method=="BFGS") {
    
    if (!is.null(hs)) {
      warning("warning: Hessian function will be ignored for quasi-Newton (non-sparse) method.")
    }
    
    if (method=="SR1") {
      if (con$preconditioner==1) {
        warning("warning:  Cannot use Cholesky decomposition as preconditioner for SR1 method.  Using identity preconditioner instead.")
        con$preconditioner <- 0L
      }
      con$quasi.newton.method <- 1L
    } else {
      con$quasi.newton.method <- 2L
    }
    ##    res <- .Call("quasiTR", x, fn1, gr1, con)    
    res <- quasiTR(x, fn1, gr1, con)
}
 
  return(res)
}

