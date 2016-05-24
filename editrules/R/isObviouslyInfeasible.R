
#' Check for obvious contradictions in a set of edits
#' 
#' Obvious contradictions are edits of the form \eqn{1 < 0}, or categorical
#' edits defining that a record fails for any value combination If
#' this function evaluates to \code{TRUE}, the set of edits is guaranteed
#' infeasible. If it evaluates to \code{FALSE} this does not garuantee feasibility.
#' See \code{\link{isFeasible}} for a complete test.
#'
#' @param E An \code{\link{editset}}, \code{\link{editmatrix}}, \code{\link{editarray}}, \code{\link[=disjunct]{editlist}} or \code{\link[=disjunct]{editenv}}
#' @param ... Arguments to be passed to or from other methods.
#' @return A \code{logical} for objects of class \code{\link{editset}}, \code{\link{editarray}} or \code{\link{editmatrix}}. 
#'  A \code{logical}  vector in the case of an \code{\link[=disjunct]{editlist}} or \code{\link[=disjunct]{editset}}.
#'
#' @export
#' @seealso \code{\link{isObviouslyRedundant}}, \code{\link{isFeasible}}
isObviouslyInfeasible <- function(E,...){
    UseMethod("isObviouslyInfeasible")
}


#'
#' @method isObviouslyInfeasible editmatrix
#' @param tol Tolerance for checking against zero.
#' @seealso \code{\link{eliminate}} \code{\link{editmatrix}}
#' @rdname isObviouslyInfeasible
#' @export
#' 
isObviouslyInfeasible.editmatrix <- function(E, tol=sqrt(.Machine$double.eps), ...){
    if ( !isNormalized(E) ) E <- normalize(E)
    A <- getAb(E)
    operators <- getOps(E)
    ib <- ncol(A)
    zeroCoef <- rowSums(abs(A[,-ib,drop=FALSE])) <= tol  
    b <- round(A[,ib],ceiling(-log10(tol)))    
    if ( any(zeroCoef & operators == "<"    &  b <= 0) || 
         any(zeroCoef & operators == "<="   &  b <  0) || 
         any(zeroCoef & operators == c("==") &  abs(b) > tol)) return(TRUE)
    return(FALSE)
}

#'
#' @method isObviouslyInfeasible editarray
#' @rdname isObviouslyInfeasible
#' @export
#' 
isObviouslyInfeasible.editarray <- function(E,...){
    any(rowSums(E)==ncol(E))
}

 
isObviouslyInfeasible.NULL <- function(E,...){
    FALSE
}


#'
#' @method isObviouslyInfeasible editset
#' @rdname isObviouslyInfeasible
#' @export
#' 
isObviouslyInfeasible.editset <- function(E,...){
    isObviouslyInfeasible(E$num) || isObviouslyInfeasible(E$mixcat)
}

#'
#' @method isObviouslyInfeasible editlist
#' @rdname isObviouslyInfeasible
#' @export
#' 
isObviouslyInfeasible.editlist <- function(E,...){
    vapply(E,isObviouslyInfeasible, FUN.VALUE=FALSE)
}

#'
#'
#'
#' @method isObviouslyInfeasible editenv
#' @rdname isObviouslyInfeasible
#' @export
#'
#'
isObviouslyInfeasible.editenv <- function(E,...){
    # note: environments are coerced to lists by lapply
    vapply(E,isObviouslyInfeasible, FUN.VALUE=FALSE)
}








