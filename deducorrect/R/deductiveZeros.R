#' Find out which variables can deductively be imputed with 0
#'
#' @param E \code{editmatrix} or Equality constraint matrix.
#' @param x named numeric vector. Naming is optional if \code{E} is an equality constraint matrix.
#' @param ... extra parameters to pass to \code{deductiveZeros,matrix}
#' @example ../examples/deductiveZeros.R
#' @export
deductiveZeros <- function(E, x,...){
    UseMethod('deductiveZeros')
}

#' Find out which variables can deductively be imputed with 0
#'
#' Interface for deductiveZeros for objects of class editmatrix. This interface
#' is robust for variables in \code{x} not occuring in \code{E}.
#'
#' @method deductiveZeros editmatrix
#' @rdname deductiveZeros
#' @export
deductiveZeros.editmatrix <- function(E, x, ...){
    eq <- getOps(E) == '=='
    vars <- getVars(E)
    xvar <- names(x)
    
    # prepare outputvector
    ddz <- logical(length(x))
    names(ddz) <- xvar
    
    # dump stuff that doesn't matter
    x <- x[xvar %in% vars]
    A <- getA(E[eq,])[ , match(xvar,vars, nomatch=0), drop=FALSE]

    # search variables s.t. nonnegativity constraints.
    nnvars <-  getVars(reduce(E[nonneg(E),]))
    nn <- logical(length(x))
    names(nn) <- names(x)
    nn[nnvars] <- TRUE

    ddz[names(x)] <- deductiveZeros.matrix(E=A,x=x, b=getb(E[eq,]), nonneg=nn, ...)
    ddz
}


#' Determine which missing values can deductively be imputed with zero
#'
#' Suppose \eqn{x} is a record under linear constraints \eqn{Ax=b} and \eqn{x\geq0}.
#' In certain cases some missing values can be imputed uniquely with zeros. For example,
#' in the case that \eqn{x_1+x_2=x_3}, if \eqn{x_2} is missing and \eqn{x_1=x_3\geq 0},
#' then \eqn{x_2} is uniquely determined to be 0. This function returns a boolean 
#' vector indicating which of the missing values are uniquely determined to be zero.
#'
#' There is some added flexibility. Users my define 'extra missings' by specifying the \code{adapt} vector.
#'
#' By default it is assumed that all values must obey the nonnegativity constraint. However this
#' can be determined by specifying \code{nonneg}. 
#'
#' @method deductiveZeros matrix
#'
#' @param b Equality constraint constant vector
#' @param adapt logical vector. Extra values to adapt, order must be the same as in \code{x}
#' @param nonneg logical vector of length(x). Determines which x-values have to obey nonnegativity constraints.
#' @param roundNearZeros Round near zero values of \code{A} before determining the sign?
#' @param tol tolerance used for zero-rounding.
#'
#' @seealso \code{\link{deduImpute}}, \code{\link{solSpace}}
#'
#' @rdname deductiveZeros
#' @export
deductiveZeros.matrix <- function(
    E, 
    x, 
    b, 
    adapt=logical(length(x)), 
    nonneg=rep(TRUE,length(x)), 
    roundNearZeros =TRUE, 
    tol=sqrt(.Machine$double.eps),
    ...
){
    m <- is.na(x) | adapt
    Amis <- E[,m,drop=FALSE]
    if (roundNearZeros) Amis[abs(Amis) < tol] <- 0
    Aobs <- E[,!m,drop=FALSE]
    bx <- b-Aobs%*%x[!m]
    r <- abs(bx) < tol & as.vector(apply(sign(Amis), 1, function(a) all(a>=0) | all(a<=0)))
    m & colSums(abs(E[r,,drop=FALSE])) > 0 & nonneg    
}



# which edits are nonnegativity constraints?
# works only on NORMALIZED editmatrix
nonneg <- function(E,tol=sqrt(.Machine$double.eps)){
    rowSums(abs(getA(E)) > tol) == 1 & abs(getb(E)) < tol & getOps(E) == '<='
}

