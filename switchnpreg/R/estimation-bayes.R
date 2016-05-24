###########################################
## Functions to estimate the parameters  ##
###########################################

fj_bayes <- function(x, y, pij, sigma2, lambda, Xdiff = outer(x, x, `-`), cov.function, ...){
    Kj <- cov.function(Xdiff, lambda = lambda, ...)
    Aj <- Kj + diag(sigma2 / pij)
    huge <- diag(Aj) > 10^8
    Ajinv <- Aj; Ajinv[,] <- 0
    Ajinv[!huge,!huge] <- solve(Aj[!huge,!huge])

    H_j <- Kj %*% Ajinv
    fjhat <- H_j %*%y
    tracej <- sum(pij*diag(H_j))
    
    list(f_hat = fjhat,
         H = H_j,
         trace = tracej)
}

## the covariance functions depends on x and on a set of parameters lambda
## in this case x will be the matrix XDiff, so that we do not need to calculate
## this all the time
covariance <- function(x, lambda, ...) {
    ## lambda is a parameter, it can be a vector depending on the
    ## covariance structure you choose;
    ## x here is a matrix with entries (x_i - x_j)
    (1/(sqrt(2*pi)*lambda))*exp(-1/2*(x^2/(lambda^2)))
}


covariance_fixed_u <- function(x, lambda, u, ...) {
    u * exp(-1/2*(x^2/(lambda^2)))
}
