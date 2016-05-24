#' Solution space for missing values under equality constraints
#'
#' @param E and \code{editmatrix} or equality constraint matrix
#' @param x a named numeric vector. 
#' @param ... Extra parameters to pass to \code{solSpace.matrix}
#' @return A \code{list} with elements \eqn{x0} and \eqn{C} or \code{NULL} if the solution space is empty
#'
#' @seealso \code{\link{deduImpute}}, \code{\link{deductiveZeros}}
#' @example ../examples/solSpace.R
#' @export
solSpace <- function(E,x,...){
    UseMethod('solSpace')
}


#' solSpace method for editmatrix
#' 
#' @param checkFeasibility Check if the observed values can lead to a consistent record
#' @method solSpace editmatrix
#' @rdname solSpace
#' @export
solSpace.editmatrix <- function(E, x, adapt=logical(length(x)), checkFeasibility=TRUE, ...){
    xvar <- names(x)
    vars <- getVars(E)
    ivar <- xvar %in% vars
    if ( checkFeasibility ){
        ii <- ivar & !(is.na(x)|adapt)
        if ( !isFeasible(substValue(E,names(x)[ii], x[ii])) ) return(NULL)
    }

    x <- x[ivar]
    eq <- getOps(E) == '=='
    if ( length(eq) == 0 ) return(NULL)
    v <- match(xvar,vars,nomatch=0)
    A <- getA(E[eq,])[,v,drop=FALSE]
    b <- getb(E[eq,])
    solSpace.matrix(A, x, b, adapt[ivar], ...)
}

#' Determine space of solutions for missing value problem.
#'
#' This function finds the space of solutions for a numerical record \eqn{x} with missing values under 
#' linear constraints \eqn{Ax=b}. Write \eqn{x=(x_{obs},x_{miss})}.
#' Then the solution space for \eqn{x_{miss}} is given by \eqn{x_0 + Cz}, where \eqn{x_0} is
#' a constant vector, \eqn{C} a constant matrix and \eqn{z} is any real vector of dimension
#' \code{ncol(C)}. This function computes \eqn{x_0} and \eqn{C}.
#'
#' The user can specify extra fields to include in \eqn{x_{miss}} by specifying \code{adapt}.
#' Also note that the method rests on the assumtion that all nonmissng values of \eqn{x} are
#' correct.
#'
#' The most timeconsuming step involves computing the generalized inverse of \eqn{A_{miss}} 
#' using \code{MASS::ginv} (code copied from MASS to avoid dependency). See the package
#' vignette and De Waal et al. (2011) for more details.
#'
#'
#' @references 
#' T. De Waal, J. Pannekoek and S. Scholtus (2011) Handbook of statistical data editing 
#' Chpt 9.2.1
#' 
#' Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with
#' S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#'
#' @method solSpace matrix
#'
#' @param b Equality constraint constant vector
#' @param adapt A named logical vector with variables in the same order as in x
#' @param tol tolerance used to determine 0-singular values when determining 
#'      generalized inverse and to round coefficients of C to zero. See \code{MASS::ginv}.
#'
#' @rdname solSpace
#' @export
solSpace.matrix <- function(
    E, 
    x, 
    b, 
    adapt=logical(length(x)),
    tol=sqrt(.Machine$double.eps),
    ...){


    m <- is.na(x) | adapt
    if ( sum(m)==0 ) return(NULL)
    Amis <- E[,m,drop=FALSE]
    Aobs <- E[,!m,drop=FALSE]
    Ainv <- ginv(Amis, tol)
    x0 <- Ainv%*%(b-Aobs%*%x[!m])
    C <-  Ainv%*%Amis - diag(1,nrow(Ainv))
    C[abs(C) < tol] <- 0
    if (!is.null(names(x))){
        rownames(x0) <-  names(x)[m]
        dimnames(C) <- list(names(x)[m],names(x)[m])
    }
    list(x0=x0, C=C)
}

#' Impute values from solution space
#'
#' Given a record \eqn{x} with observerd \eqn{x_{obs}} and missing values \eqn{x_{miss}} under
#' linear equality constraints \eqn{Ax=b}. The function \code{\link{solSpace}} returns
#' the solution space which can be written as \eqn{x_{miss} = x_0 + Cz}, where \eqn{x_0} is
#' are a constant vector (of dimension d=\code{length}\eqn{(x_{miss})}) and \eqn{C} a constant
#' matrix of dimension \eqn{d\times d}. 
#'
#' If \eqn{C} has rows equal to zero, then those missing values may be imputed deductively.
#' For the other missing values, some \eqn{z} must be chosen or another imputation method
#' used.
#'
#' The function \code{imputess} imputes missing values in a vector \eqn{x}, based on the
#' solution space and some chosen vector \eqn{z}. If no \eqn{z} is passed as argument, only
#' deductive imputations are performend (i.e. some missings may be left).
#' 
#' If \eqn{C} is a named matrix (as returned by \code{\link{solSpace}}), rows of \eqn{x0} and \eqn{C}
#' are matched by name to \eqn{x}. Otherwise it is assumed that the missings in \eqn{x} occur in the order
#' of the rows in \eqn{C} (which is also the case when x0 and C are computed by \code{\link{solSpace}}).
#'
#' @param x (named) numerical vector to be imputed
#' @param x0 \code{x0} outcome of \code{\link{solSpace}}
#' @param C \code{C} outcome of \code{\link{solSpace}}
#' @param z real vector of dimension \code{ncol(C)}.
#' @param tol tolerance used to check which rows of \code{C} equal zero.
#'
#' @example ../examples/deduImpute.R
#'
#' @export
imputess <- function(x, x0, C, z=NULL, tol=sqrt(.Machine$double.eps)){
    imiss <- colnames(C)
    if ( is.null(imiss) ) imiss <- is.na(x)
    if (is.null(z)){ # just impute unique values
        u <- rowSums(abs(C)) < tol
        x[imiss][u] <- x0[u]
    } else {
        x[imiss] <- x0 +  C%*%z
    }
    x
}



# Generalized matrix inverse. Code copied from MASS:::ginv 
#
#
ginv <- function (X, tol = sqrt(.Machine$double.eps)) {
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}



