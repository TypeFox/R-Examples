rpdatasim <- function(n, replications, rdist, ...) {
    stopifnot(is.function(rdist), is.numeric(n), is.numeric(replications))

    matrix(data = rdist(n * replications, ...),
           nrow = n,
           ncol = replications)
}
