#' Find obvious redundancies in set of edits
#' 
#' Detect simple redundancies such as duplicates or edits of the form \code{0 < 1} or \code{0 == 0}.
#' For categorical edits, simple redundancies are edits that define an empty subregion
#' of the space of all possible records (no record can ever be contained in such a region).
#'
#' @param E An \code{\link{editset}}, \code{\link{editmatrix}}, \code{\link{editarray}}, 
#'      \code{\link[=disjunct]{editlist}} or \code{\link[=disjunct]{editenv}}
#' @param duplicates \code{logical}: check for duplicate edits? For an \code{\link{editset}}, 
#'      \code{\link[=duplicated]{editlist}} or \code{\link[=duplicated]{editenv}} this should be a logical
#'      2-vector indicating which of the numerical or categorical edits should be checked for duplicates.
#' @param ... parameters to be passed to or from other methods. 
#' 
#' @return logical vector indicating which edits are (obviously) redundant
#' @seealso \code{\link{isObviouslyInfeasible}}, \code{\link{isSubset}}
#' @export 
isObviouslyRedundant <- function(E, duplicates=TRUE, ...){
    UseMethod("isObviouslyRedundant")
}


# @method isObviouslyRedundant matrix
#
# @param operators character vecor of comparison operators in \code{<, <=, ==} of length \code{nrow(E)}
# @param tol tolerance to check for zeros.
#
# @rdname isObviouslyRedundant
# @keywords internal
# @seealso \code{\link{isObviouslyRedundant}}, \code{\link{isObviouslyRedundant.editmatrix}}
isObviouslyRedundant.matrix <- function(
    E, 
    operators, 
    tol=sqrt(.Machine$double.eps), 
    ... ){
    ib <- ncol(E)
    zeroCoef <- rowSums(abs(E[,-ib,drop=FALSE])) < tol
    as.vector(
        zeroCoef & ( (operators %in% c("==","<=")  & abs(E[,ib]) < tol) 
                   | (operators %in% c("<", "<=")  & E[,ib] > tol)
                   )
    )
}


#' @method isObviouslyRedundant editmatrix
#' @rdname isObviouslyRedundant
#' 
#' @export
isObviouslyRedundant.editmatrix <- function(E, duplicates=TRUE, ...){
    if ( !isNormalized(E) ) E <- normalize(E)
    I <- isObviouslyRedundant.matrix(getAb(E), operators=getOps(E), ...)
    if ( duplicates ) I <- I | duplicated.editmatrix(E)
    I
}


#' @method isObviouslyRedundant editarray
#' @rdname isObviouslyRedundant
#' @export
isObviouslyRedundant.editarray <- function(E, duplicates=TRUE, ...){
    if ( ncol(E) == 0 ) return(logical(0))
    if ( ncol(E) == 1 ) return(as.vector(E))
    ind <- getInd(E)
    red <- isRedundant.boolmat(getArr(E),getInd(E))
    if ( duplicates ) red <- red | duplicated.editarray(E)
    red
}


# Check redundancy in editarray after disection
#
# @keywords internal
isRedundant.boolmat <- function(A, ind){
    if ( nrow(A) == 1 ) return(any(vapply(ind,function(i) sum(A[,i])==0,FUN.VALUE=TRUE)))
    apply(
        vapply(ind, function(i) rowSums(A[,i,drop=FALSE])==0, FUN.VALUE=logical(nrow(A))),
        1,any
    )
}

#
# For an \code{\link{editset}} it checks for obvious redundancies in the numerical
# and categorical/mixed parts separately. Arguments \code{duplicates} must be a
# logical 2-vector, the first element corresponding to the numeric part, the second 
# element to the conditional part.
#

#' @method isObviouslyRedundant editset
#' @rdname isObviouslyRedundant
#' @export
isObviouslyRedundant.editset <- function(E, duplicates=rep(TRUE,2), ...){
    c(
        isObviouslyRedundant(E$num, duplicates=duplicates[1], ...),
        isObviouslyRedundant(E$mixcat, duplicates=duplicates[2], ...)
    ) 
}



#
# 
# For an \code{\link[=disjunct]{editlist}} or \code{\link[=disjunct]{editenv}},
# a list of logical vectors is returned.
#
#' @method isObviouslyRedundant editlist
#' @rdname isObviouslyRedundant
#' @export
isObviouslyRedundant.editlist <- function(E, duplicates=rep(TRUE,2), ...){
    lapply(E, isObviouslyRedundant.editset, duplicates=duplicates, ...)
}


# 
# For an \code{\link[=disjunct]{editlist}} or \code{\link[=disjunct]{editenv}},
# a list of logical vectors is returned.
#

#' @method isObviouslyRedundant editenv
#' @rdname isObviouslyRedundant
#' @export
isObviouslyRedundant.editenv <- function(E, duplicates=rep(TRUE,2), ...){
    lapply(E, isObviouslyRedundant.editset, duplicates=duplicates, ...)
}



