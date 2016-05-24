multmix.init <- function(y, lambda = NULL, theta = NULL, k = 2){

n <- nrow(y)
p <- ncol(y)

    if (is.null(theta)) {
        theta = matrix(runif(p * k), k, p)
        theta = theta/apply(theta, 1, sum)
    }
    else k = nrow(theta)
    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    } else k = length(lambda)

list(lambda=lambda, theta=theta, k=k)


}