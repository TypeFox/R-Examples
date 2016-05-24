
#' Adjust a data to meet linear (in)equality constraints
#'
#' Adjust a vector \eqn{\boldsymbol{x}} to meet
#' constraints \eqn{\boldsymbol{Ax} \leq \boldsymbol{b}}. 
#' 
#'
#' @param object an \code{R} object describing constraints (see details)
#' @param ... Arguments to be passed to other methods
#'
#' @return Object of class \code{\link{adjusted}}.
#'
#' @section Details:
#' \code{adjust} is a generic function allowing several definitions of the constraints in \code{object}.
#'
#' \itemize{
#'  \item[editmatrix]{If \code{object} is an \code{editmatrix}, the function will 
#'    try to match the names of \code{x} to the variable names in \code{object} before further
#'    processing. In that case the \code{length} of \code{x} is unimportant, as long as all variables in \code{object} 
#'    are also in \code{x}. Depending on the choice of \code{method}, \code{object} is converted to \code{matrix} or 
#'    \code{sparseConstraints} format before solving the adjustment problem.
#'  }
#'  \item[matrix]{If \code{object} is a \code{matrix}, you also need to provide the constant vector
#'   \code{b} and the number of equations \code{neq} to define the problem. It is assumed that the first
#'   \code{neq} rows of \code{object} and the first \code{new} elements of \code{b} correspond to equalities. No names are matched, so \code{x}
#'    must be in the correct order and must be of the right dimension.
#'    See \code{\link{sparseConstraints}} on how to translate
#'    a \code{matrix} problem to the sparse version.
#'  } 
#'
#' \item[sparseConstraints] {If \code{object} is of class \code{\link{sparseConstraints}}, 
#'   the sparse method is used to adjust \code{x}. Some basic checks on \code{x} and \code{w} 
#'   are performed, but no attempt is made to match names of \code{x} to those of \code{object}.
#' }
#'}
#'
#' The tolerance \code{tol} is defined as the maximum absolute value of the difference vector 
#' \eqn{\boldsymbol{Ax}-\boldsymbol{b}} for equalities. For inequalities, the difference vector
#' is set to zero when it's value is lesser than zero (i.e. when the restriction is obeyed). The
#' function keeps iterating until either the tolerance is met, the number of allowed iterations is
#' exceeded or divergence is detected. 
#'
#' @section Note:
#' \code{adjust} does not perform any consistency checks. When the system of constraints is
#' contradictory (\emph{e.g.} \eqn{x>1} and \eqn{x<0}) this will result in either divergence
#' or in exceeding the number of iterations. 
#'
#'
#'
#' @example ../examples/adjust.R
#' @export
adjust <- function(object, ...){
   UseMethod('adjust')
}

#'
#' @method adjust editmatrix
#' @param method use dense or sparse matrix method.
#' @export
#' @rdname adjust
adjust.editmatrix <- function(object, x, w=rep(1,length(x)), method=c('dense','sparse'), ...){
  method <- match.arg(method)
   if (!isNormalized(object)) object <- normalize(object)
	object <- reduce(object)
   # match names 
   if ( !is.null(names(x)) ){
      J <- match(getVars(object), names(x))
   } else {
      stopifnot(length(x) == length(getVars(object)))
      J <- 1:length(x)
   }
   u <- x[J]
   w <- w[J]
   ops <- getOps(object)
   I <- order(ops,decreasing=TRUE)
   neq <- sum(ops == "==")


	if ( method == 'sparse' ){
		y <- adjust.sparseConstraints(
         sparseConstraints(object), 
         x = u, 
         w = w, 
         ...
      )
	} else {
		y <- adjust.matrix(
			object = getA(object)[I,,drop=FALSE], 
			b      = getb(object)[I], 
			x      = u,
			neq    = neq,
         w      = w,
			... 
		)
	} 
		
   
   x[J] <- y$x
   y$x <- x
   y
}

#'
#' @method adjust sparseConstraints
#' @export
#' @rdname adjust
adjust.sparseConstraints <- function(object, x, w=rep(1.0,length(x)), tol=1e-2, maxiter=1000L, ...){

   stopifnot(
		is.numeric(x),
		length(x) == object$.nvar(),
      all_finite(x),
      length(w) == length(x),
      is.numeric(w),
      all_finite(w),
      all(w>0),
      tol > 0,
      is.finite(tol),
      maxiter > 0,
      is.finite(maxiter)
   )

   y <- object$.adjust(x, w, tol, maxiter)
   
}



#' 
#' @param b Constant vector of the constraint system \eqn{Ax\leq b}
#' @param x The vector to be adjusted
#' @param neq the first \code{neq} linear relations are equalities.
#' @param w A positive weight vector
#' @param tol The maximum allowed deviation from the constraints (see details).
#' @param maxiter maximum number of iterations
#'
#' @method adjust matrix
#' @export
#' @rdname adjust
adjust.matrix <- function(object, b, x, neq=length(b), w=rep(1.0,length(x)), tol=1e-2, maxiter=1000L, ...){
   stopifnot(
		is.numeric(x),
		length(x) == ncol(object),
      all_finite(x),
		is.numeric(b),
		length(b) == nrow(object),
      all_finite(b),
      length(w) == length(x),
      is.numeric(w),
      all_finite(w),
      all(w>0),
      tol > 0,
      is.finite(tol),
      maxiter > 0,
      is.finite(maxiter)
   )


   storage.mode(object) <- "double"
   
   
   t0 <- proc.time()
   y <- .Call("R_dc_solve", 
      object, 
      as.double(b), 
      as.double(w),
      as.integer(neq),
      as.double(tol),
      as.integer(maxiter),
      as.double(x)
   )
   t1 <- proc.time()
   objective <- sqrt(sum(w*(x-as.vector(y))^2))
   new_adjusted(y, t1-t0,"dense", objective, colnames(object))
} 







