#' Generate sparse set of constraints.
#'
#' @param x R object to be translated to sparseConstraints format.
#' @param ... options to be passed to other methods
#'
#' @return Object of class \code{sparseConstraints} (see details).
#' 
#' @section Details:
#' 
#' The \code{sparseConstraints} objects holds the system \eqn{\boldsymbol{Ax}\leq \boldsymbol{b}}
#' in column sparse format, outside of \code{R}'s memory. In \code{R}, it is a \emph{reference object}.
#' In particular, it is meaningless to
#' \itemize{
#'    \item{Copy the object. You only will only generate a pointer to physically the same object.}
#'    \item{Save the object. The physical object is destroyed when \code{R} closes, or when \code{R}'s
#'      garbage collector cleans up a removed \code{sparseConstraints} object.}
#' }
#'
#' @export
#' @example ../examples/sparseConstraints.R
sparseConstraints = function(x, ...){
    UseMethod("sparseConstraints")

}


#'
#' @method sparseConstraints editmatrix
#' @param tol Tolerance for testing where coefficients are zero
#' @rdname sparseConstraints
#' @export
sparseConstraints.editmatrix = function(x, tol=1e-8, ...){
   if (!isNormalized(x)) normalize(x)
   x <- reduce(x,tol=tol)
   ieq <- getOps(x) == '=='
   I <- c(which(ieq),which(!ieq))
   x <- x[I,];
   e <- new.env();
   A <- getA(x);
   storage.mode(A) <- "double"
   e$.sc <- .Call("R_sc_from_matrix", A, as.double(getb(x)), as.integer(sum(ieq)), as.double(tol))
   e$.vars <- getVars(x)
   make_sc(e)
}

#'
#'
#' @method sparseConstraints matrix
#' @rdname sparseConstraints
#' @export
sparseConstraints.matrix <- function(x, b, neq=length(b), tol=1e-8,...){
	stopifnot(
		all_finite(x),
		is.numeric(b),
		all_finite(b),
		length(b) == nrow(x),
		is.numeric(neq),
		is.finite(neq),
		neq > 0,
		neq <= length(b),
		is.numeric(tol),
		is.finite(tol),
		tol > 0
	)
	storage.mode(x) <- "double"
	e <- new.env()
   e$.sc <- .Call("R_sc_from_matrix", x, as.double(x), as.integer(neq), as.double(tol))
   e$.vars <- colnames(x)
	make_sc(e)
	
}



#'
#' @method sparseConstraints data.frame
#'
#' @param b Constant vector
#' @param neq The first \code{new} equations are interpreted as equality constraints, the rest as '<='
#' @param base are the indices in \code{x[,1:2]} base 0 or base 1?
#' @param sorted is \code{x} sorted by the  first column?
#' @export
#' @rdname sparseConstraints
sparseConstraints.data.frame <- function(x, b, neq=length(b), base=min(x[,2]), sorted=FALSE, ...){
   if (length(b) != length(unique(x[,1]))){
      stop("length of b unequal to number of constraints")
   }
	
	stopifnot(
		is.numeric(x[,1]),
		all_finite(x[,1]),
		is.numeric(x[,2]),
		all_finite(x[,2]),
      all(x[,2]>=base),
		is.numeric(b),
		all_finite(b),
      is.numeric(neq),
		is.finite(neq),
		neq <= length(b),
		base %in% c(0,1)
	)


	if ( !sorted ) x <- x[order(x[,1]),,drop=FALSE]
   e <- new.env()
   e$.sc <- .Call("R_sc_from_sparse_matrix", 
      as.integer(x[,1]), 
      as.integer(x[,2]-base),
      as.double(x[,3]), 
      as.double(b),
      as.integer(neq)
   )
   make_sc(e)

}




#' 
#' @method print sparseConstraints
#' @param range integer vector stating which constraints to print
#'
#' @export
#' @rdname sparseConstraints
print.sparseConstraints <- function(x, range=1L:10L, ...){
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
   e$.adjust <- function(x, w, tol, maxiter){
      t0 <- proc.time() 
      y <- .Call('R_solve_sc_spa',
         e$.sc, 
         as.double(x), 
         as.double(w), 
         as.double(tol), 
         as.integer(maxiter)
      )
      t1 <- proc.time()
      objective <- sqrt(sum((x-as.vector(y))^2*w))
      new_adjusted(y,t1-t0, "sparse", objective, e$.vars)
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

   
   structure(e,class="sparseConstraints")
}






