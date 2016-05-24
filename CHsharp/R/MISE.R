MISE <- function(x, xgrid, sigma2, lambda, h, g, A, B) {
    II <- diag(rep(1, length(x)))
    inverse <- solve(II+lambda*B%*%t(B))
    gvector <- g(x)
    if (length(sigma2) == 1) {
        S <- diag(rep(sigma2, length(x)))
    } else {
        S <- diag(sigma2)
    } 
    variance <- diag(t(A)%*%inverse%*%S%*%inverse%*%A)
    Astar <- t(A)%*%inverse
    bias <- g(xgrid) - Astar%*%gvector
    c(var=mean(variance)*diff(range(xgrid)),
bias2 = diff(range(xgrid))*mean(bias^2), MISE = diff(range(xgrid))*
(mean(variance) + mean(bias^2)))
}


