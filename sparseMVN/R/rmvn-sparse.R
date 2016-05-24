## Copyright (C) 2013-2015 Michael Braun
#' @title Multivariate normal functions with sparse covariance/precision matrix.
#' @aliases dmvn.sparse rmvn.sparse
#' @description Efficient sampling and density calculation from a multivariate
#' normal,
#' when the covariance or precision matrix is sparse. These functions are
#' designed for MVN samples of very large dimension.
#'
#' @param n number of samples
#' @param x numeric matrix, where each row is an MVN sample.
#' @param mu mean (numeric vector)
#' @param CH An object of class dCHMsimpl or dCHMsuper that represents
#' the Cholesky factorization of either the precision (default) or covariance
#' matrix.  See details.
#' @param prec If TRUE, CH is the Cholesky decomposition of the precision
#' matrix.  If false, it is the decomposition for the covariance matrix.
#'
#' @section Details:
#' These functions uses sparse matrix operations to sample from, or compute the
#' log density of, a multivariate normal distribution  The user must compute
#' the Cholesky decomposition first, using the Cholesky function in the Matrix
#' package.  This function operates on a sparse symmetric matrix, and returns
#' an object of class dCHMsimpl or dCHMsuper (this depends on the algorithm
#' that was used for the decomposition).  This object contains information about
#' any fill-reducing permutations that were used to preserve sparsity.  The
#' rmvn.sparse and dmvn.sparse functions use this permutation information, even
#' if pivoting was turned off.
#'
#' @examples
#'    require(Matrix)
#'    m <- 20
#'    p <- 2
#'    k <- 4
#' 
#' ## build sample sparse covariance matrix
#'    Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(m)))
#'    Q2 <- cBind(Q1,Matrix(0,m*p,k))
#'    Q3 <- rBind(Q2,cBind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
#'    V <- tcrossprod(Q3)
#'    CH <- Cholesky(V)
#' 
#'    x <- rmvn.sparse(10,rep(0,p*m+k),CH, FALSE)
#'  ##  print(x)
#' 
#'    y <- dmvn.sparse(x[1,],rep(0,p*m+k), CH, FALSE)
#'  ##  print(y)
#'
#' @export
rmvn.sparse <- function(n, mu, CH, prec=TRUE) {

    if (is.na(match(class(CH),c("dCHMsimpl","dCHMsuper")))) {
        stop("CH must be an object of class 'dCHMsimpl' or 'dCHMsuper'")
    }

    k <- length(mu)

    if (!(k>0)) {
        stop("mu must have positive length")
    }
   
    if (!(n>0)) {
        stop("n must be positive")
    }

    if (!(k==dim(CH)[1])) {
        stop("dimensions of mu and CH do not conform")
    }

    if (!is.logical(prec)) {
        stop("prec must be either TRUE or FALSE")
    }

    x <- rnorm(n*k)
    dim(x) <- c(k,n)

    A <- expand(CH)

    if (prec) {
        y <- solve(Matrix::t(A$L),x) ## L'y = x
    } else {
        y <- A$L %*% x
    }

    y <- as(Matrix::crossprod(A$P,y),"matrix") ## P' %*% y

    y <- y + mu

    return(t(y))

}

#' @rdname rmvn.sparse
#' @export
dmvn.sparse <- function(x, mu, CH, prec=TRUE) {


    if (is.vector(x) | (is.atomic(x) & NCOL(x)==1)) {
        x <- matrix(x,nrow=1)
    }

    k <- length(mu)
    n <- NROW(x)

    if (!(k>0)) {
        stop("mu must have positive length")
    }

    if (!(k==dim(CH)[1])) {
        stop("dimensions of mu and CH do not conform")
    }

    if (k!=NCOL(x)) {
        stop("x must have same number of columns as the length of mu")
    }


    if (!is.logical(prec)) {
        stop("prec must be either TRUE or FALSE")
    }

    A <- expand(CH)

    detL <- sum(log(Matrix::diag(A$L)))
    C <- -0.918938533204672669541*k ## -k*log(2*pi)/2


    xmu <- t(x)-mu

    z <- as.matrix(A$P %*% xmu)

    if (prec) {
        y <- Matrix::crossprod(A$L,z)  ## L' %*% x
        log.dens <- C + detL - Matrix::colSums(y*y)/2
    } else {
        y <- solve(A$L, z) ## Ly = x
        log.dens <- C - detL - Matrix::colSums(y*y)/2
    }

    return(as.numeric(log.dens))
}

