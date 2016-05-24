#' Generate sparse set of constraints.
#'
#' Generate a constraint set to be used by \code{\link{sparse_project}}
#'
#'
#' @param object R object to be translated to sparse_constraints format.
#' @param ... options to be passed to other methods
#'
#' @return Object of class \code{sparse_constraints} (see details).
#' 
#' @section Note:
#' 
#' As of version 0.1.1.0, \code{sparseConstraints} is deprecated. Use \code{sparse_constraints} instead.
#'
#' @section Details:
#' 
#' The \code{sparse_constraints} objects holds coefficients of
#' \eqn{\boldsymbol{A}} and \eqn{\boldsymbol{b}} of the system
#' \eqn{\boldsymbol{Ax}\leq \boldsymbol{b}} in sparse format, outside of
#' \code{R}'s memory. It can be reused to find solutions for vectors to adjust.
#'
#' In \code{R}, it is a \emph{reference object}. In particular, it is meaningless to
#' \itemize{
#'    \item{Copy the object. You only will only generate a pointer to physically the same object.}
#'    \item{Save the object. The physical object is destroyed when \code{R} closes, or when \code{R}'s
#'      garbage collector cleans up a removed \code{sparse_constraints} object.}
#' }
#'
#' @section The \code{$project} method:
#' 
#' Once a \code{sparse_constraints} object \code{sc} is created, you can reuse it to optimize
#' several vectors by calling \code{sc$project()} with the following parameters:
#' \itemize{
#'   \item{\code{x}: \code{[numeric]} the vector to be optimized}
#'   \item{\code{w}: \code{[numeric]} the weight vector (of \code{length(x)}). By default all weights equal 1.}
#'   \item{\code{eps}: \code{[numeric]} desired tolerance. By default \eqn{10^{-2}} }
#'   \item{\code{maxiter}: \code{[integer]} maximum number of iterations. By default 1000.}
#' }
#' The return value of \code{$spa} is the same as that of \code{\link{sparse_project}}.
#' 
#' @seealso \code{\link{sparse_project}}, \code{\link{project}}
#' @export
#' @example ../examples/sparse_constraints.R
sparse_constraints <- function(object, ...){
  UseMethod("sparse_constraints")
}

#' @rdname sparse_constraints
#' @export
sparseConstraints <- function(object,...){
  .Deprecated("sparse_constraints")
  sparse_constraints(object,...)
}


#' Read sparse constraints from a \code{data.frame}
#' 
#' @method sparse_constraints data.frame
#'
#' @param b Constant vector
#' @param neq The first \code{new} equations are interpreted as equality constraints, the rest as '<='
#' @param base are the indices in \code{object[,1:2]} base 0 or base 1?
#' @param sorted is \code{object} sorted by the  first column?
#' @export
#' @rdname sparse_constraints
sparse_constraints.data.frame <- function(object, b, neq=length(b), base=1L, sorted=FALSE, ...){

  if (length(b) != length(unique(object[,1]))){
    stop("length of b unequal to number of constraints")
  }
	
  stopifnot(
    is.numeric(object[,1])
    , all_finite(object[,1])
    , is.numeric(object[,2])
    , all_finite(object[,2])
    , all(object[,2]>=base)
    , is.numeric(b)
    , all_finite(b)
    , is.numeric(neq)
    , is.finite(neq)
    , neq <= length(b)
    , base %in% c(0,1)
  )


	if ( !sorted ) object <- object[order(object[,1]),,drop=FALSE]
   e <- new.env()
   e$.sc <- .Call("R_sc_from_sparse_matrix", 
      as.integer(object[,1]), 
      as.integer(object[,2]-base),
      as.double(object[,3]), 
      as.double(b),
      as.integer(neq)
   )
   make_sc(e)

}




#' Print sparse_constraints object
#' 
#' @method print sparse_constraints
#' @param range integer vector stating which constraints to print
#' @param x an object of class \code{sparse_constraints}
#' @export
#' @rdname sparse_constraints
print.sparse_constraints <- function(x, range=1L:10L, ...){
   x$.print()
}

# e: environment containing an R_ExternalPtr
make_sc <- function(e){
   #
   
  e$.pointer <- function(){
    e$.sc
  }
  
  e$.nvar <- function(){
    .Call("R_get_nvar", e$.sc)
  }
  
  e$.nconstr <- function(){
    .Call("R_get_nconstraints", e$.sc)
  }
  
  e$.print <- function(range){
    if ( missing(range) & e$.nvar() > 10 ) range = numeric(0)
    if ( missing(range) & e$.nvar() <=10 ) range = 1L:10L
    vars = e$.vars
    if ( is.null(vars) ) vars = character(0);
  
    stopifnot(all(range >= 1))
    range = range-1;
  
    dump <- .Call("R_print_sc",e$.sc, vars, as.integer(range))
  }

  # adjust input vector minimally to meet restrictions.
  e$project <- function(x, w=rep(1,length(x)), eps=1e-2, maxiter=1000L){
    stopifnot(
      eps > 0
      , maxiter > 0
      , all_finite(w)
      , all_finite(x)
    )
    t0 <- proc.time() 
    y <- .Call('R_solve_sc_spa',
       e$.sc, 
       as.double(x), 
       as.double(w), 
       as.double(eps), 
       as.integer(maxiter)
    )
    t1 <- proc.time()
    objective <- sqrt(sum((x-as.vector(y))^2*w))
    
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

  e$.diffsum <- function(x){
    stopifnot(length(x)==e$.nvar())
    .Call("R_sc_diffsum", e$.sc, as.double(x)) 
  }
  
  e$.diffmax <- function(x){
    stopifnot(length(x)==e$.nvar())
    .Call("R_sc_diffmax", e$.sc, as.double(x)) 
  }
  
  e$.multiply <- function(x){
    stopifnot(length(x) == e$.nvar());
    .Call("R_sc_multvec", e$.sc, as.double(x))
  }
  
  e$.diffvec <- function(x){
    stopifnot(length(x) == e$.nvar())
    .Call("R_sc_diffvec", e$.sc, as.double(x))
  }
  
  
  structure(e,class="sparse_constraints")
}






