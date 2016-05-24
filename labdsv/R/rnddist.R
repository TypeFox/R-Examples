rnddist <- function (size, method='metric', sat = 1.0, upper=FALSE, diag=FALSE)
{
    tmp <- matrix(runif(size^2),ncol=size)
    x <- as.dist(tmp)
    if (method == 'metric') {
        y <- metrify(x, upper=upper, diag=diag)
        attr(y, "method") <- "metric random"
    } else if (method == 'euclidean') {
        y <- euclidify(x, upper=upper, diag=diag)
        attr(y,"method") <- 'euclidean random'
    } else {
        stop("you must specify 'metric' or 'euclidean' as the method")
    }

    y <- y/max(y)
        if (sat != 1.0) {
        y <- y / sat
        y[y>1.0] <- 1.0
    }
    attr(y, "call") <- match.call()
    attr(y, "Diag") <- diag
    attr(y, "Upper") <- upper
    y
}
