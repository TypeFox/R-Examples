### Returns the lambda matrix, $\Lambda$, to calculate the first
### divided differences simplicity measure of a vector $v \in R^d$ via
### $v^T \Lambda\ v$.
###
### Arguments:
###   x - vector containing ordered values of the environment levels
lambda_first <- function(x) {
    d <- length(x)
    
    D <- mdiff(d)
    W <- weight_first(x)
    
    4 * diag(d) - crossprod(W %*% D)
}


### Returns the lambda matrix, $\Lambda$, to calculate the second
### divided differences simplicity measure of a vector $v \in R^d$ via
### $v^T \Lambda\ v$.
###
### Arguments:
###   x - vector containing ordered values of the environment levels
lambda_second <- function(x) {
    d <- length(x)
    
    D1 <- mdiff(d)
    D2 <- mdiff(d-1)
    
    W1 <- weight_second1(x)
    W2 <- weight_second2(x)
    
    lambda <- crossprod(W2 %*% D2 %*% W1 %*% D1)
    max(eigen(lambda, only.values = TRUE)$values) * diag(d) - lambda
}


mdiff <- function(d) {
    x <- matrix(0, d-1, d)
    diag(x) <- -1
    diag(x[ , seq(2, d)]) <- 1
    x
}


weight_first <- function(x) {
    w <- diff(x)
    minw <- min(w)
    diag(sqrt(minw/w))
}


weight_second1 <- function(x) {
    w <- diff(x)
    diag(1/w)
}


weight_second2 <- function(x) {
    w <- diff(x, lag=2)
    diag(sqrt(1/w))
}
