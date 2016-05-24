.isSingular <- function(mat)
{
##    return( - (determinant(mat, logarithm = TRUE)$modulus[1] - 0)/ncol(mat) > 50)
    p <- ncol(mat)
    if(!is.qr(mat))
        mat <- qr(mat)
    return(mat$rank < p)
}

.check_vars_numeric <- function(mf)
{
    ## we need to test just the columns which are actually used.
    mt <- attr(mf, "terms")
    mterms <- attr(mt, "factors")
    mterms <- rownames(mterms)[apply(mterms, 1, any)]
    any(sapply(mterms, function(x) is.factor(mf[,x]) || !is.numeric(mf[,x])))
}

vecnorm <- function(x, p=2) sum(x^p)^(1/p)

##
## Several Matlab-like utility functions ======================================================
##
##
## Return the square root of a symetric positive definite matrix
sqrtm <- function(A){
    ##
    ## [E D] = eig(A); sqrtm(A) = E * sqrt(D) * E'
    ##
    if(!is.matrix(A) || ncol(A) != nrow(A))
        stop("The matrix A must be a square matrix\n")

    ee <- eigen(A)
    if(any(ee$values < 0)) {
        stop("The matrix A must be positive definite.")
    }

    ee$vectors %*% diag(sqrt(ee$values)) %*% t(ee$vectors)
}

## Return an n by p matrix of ones
##ones <- function(n=1, p=1){
##    matrix(1, nrow=n, ncol=p)
##}

## Return an n by p matrix of zeros
##zeros <- function(n=1, p=1){
##    matrix(0, nrow=n, ncol=p)
##}

##  <Matlab>
##      a=[1 2 ; 3 4];
##      repmat(a,2,3)
##
##  <R>
##      a <- matrix(1:4,2,byrow=T)
##      repmat(a,2,3)
##
##
##  a <- 1:4; n=10
##  matrix(rep(a, times=n), nrow=n, byrow=TRUE)
##

##repmat <- function(A, n, p) {
##
##    if(is.vector(A))    # we need a column matrix, not a vector, speaking in R terms
##        A <- t(A)
##    kronecker(matrix(1,n,p), A)
##}
