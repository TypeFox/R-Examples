## Gaussian process regression for bdpopt

## Construct a specific entry of the squared exponential covariance matrix
squared.exp.cov <- function(x1, x2, sigma.f, lengths) {   
    sigma.f^2 * exp( -(1 / 2) * sum( ( (x1 - x2) / lengths )^2 ) )
}

## Construct a squared exponential covariance matrix, using the fact that it is symmetric
squared.exp.cov.matrix <- function(X, sigma.f, lengths) {
    n <- ncol(X)
    A <- matrix(vector(mode = "numeric"), nrow = n, ncol = n)

    for (i in 1:n) {
        A[i,i:n] <- apply(X[,i:n, drop = FALSE], 2,
                          function(x) sum( ( (X[,i] - x) / lengths )^2 ) )        
    }

    A[lower.tri(A)] = t(A)[lower.tri(A)]
    sigma.f^2 * exp( -(1 / 2) * A )   
}

## Compute the Cholesky factorisation, stopping if the covariance matrix is not positive definite
cholesky.factorisation <- function(K, vars) {    
    U <- try(chol(K + diag(vars)))
    if (inherits(U, "try-error"))
        stop("Cholesky factorisation failed for GPR regression, you might try loess regression instead")
    
    U
}

## Gaussian process regression using a squared exponential covariance function.
## Given a matrix of training points X (each column a point) and a vector y of responses,
## construct and return a function giving the posterior predictive mean corresponding to a new observation.
## Format for hyperparams: c(sigma.f, length.1, ..., length.d),
## d being the dimension of the choice space for the original problem.
## vars is a vector of variances for the observational errors, which are assumed to be normally distributed.
gpr.given.hyperparams <- function(X, y, vars, sigma.f, lengths) {   
    ## Find the upper triangular part of the Cholesky factorisation
    K <- squared.exp.cov.matrix(X, sigma.f, lengths)
    U <- cholesky.factorisation(K, vars)
    alpha <- backsolve(U, forwardsolve(t(U), y))
    
    X.cols <- mat.cols(X)
    k.star <- function(x) vapply(X.cols, function(x2) squared.exp.cov(x, x2, sigma.f, lengths), 0)

    function(x) sum(k.star(x) * alpha)    
}

## The objective function to minimize, corresponding to maximisation of the logarithm
## of the marginal likelihood for a squared exponential model.
squared.exp.lh.obj <- function(X, y, vars, sigma.f, lengths) {             
    ## Find the upper triangular part of the Cholesky factorisation
    K <- squared.exp.cov.matrix(X, sigma.f, lengths)
    U <- cholesky.factorisation(K, vars)    
    alpha <- backsolve(U, forwardsolve(t(U), y))    
    
    sum(y * alpha + 2 * log(diag(U)))
}

## The gradient of the objective function above.

## Since k(x1, x2) = sigma.f^2 * exp( - (1 / 2) * sum( ( (x1 - x2) / lengths )^2 ) ),
## letting lengths = (l.i), i = 1, ..., d,

## d k(x1, x2) / d sigma.f = 2 sigma.f * exp( - (1 / 2) * sum( ( (x1 - x2) / lengths )^2 ) )
## ---> d K / d sigma.f = 2 K / sigma.f

## d k(x1, x2) / d l.i =
## sigma.f^2 * exp( - (1 / 2) * sum( ( (x1 - x2) / lengths )^2 ) ) * (- (1 / 2) * (-2) * (x1[i] - x2[i])^2 / l.i^3) =
## k(x1, x2) * (x1[i] - x2[i])^2 / l.i^3   

squared.exp.lh.obj.grad <- function(X, y, vars, sigma.f, lengths) {
    ## Find the upper triangular part of the Cholesky factorisation
    K <- squared.exp.cov.matrix(X, sigma.f, lengths)
    U <- cholesky.factorisation(K, vars)    
    
    alpha <- backsolve(U, forwardsolve(t(U), y))
    K.inv <- chol2inv(U)
    K.alpha <- K.inv - alpha %*% t(alpha)
    
    ## Compute the partial derivative w.r.t. sigma.f        
    d.K.d.sigma.f <- 2 * K / sigma.f
    d.f.d.sigma.f <- tr(K.alpha %*% d.K.d.sigma.f)      
    
    ## Compute the partial derivatives w.r.t. the lengths            
    A.list <- lapply(mat.rows(X), function(row) do.call(cbind, lapply(row, function(x) (row - x)^2)))
    
    d.K.d.lengths <- mapply(function(A, l) tr(K.alpha %*% (K * A / l^3)), A.list, lengths)                
    
    ## Return the gradient
    c(d.f.d.sigma.f, d.K.d.lengths)
}
