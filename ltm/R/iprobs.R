iprobs <-
function (betas, z) {
    n <- length(z)
    gammas <- lapply(betas, function (x) {
        nx <- length(x)
        cbind(plogis(matrix(x[-nx], n, nx - 1, TRUE) - x[nx] * z), 1)
    })
    lapply(gammas, function (x) {
        nc <- ncol(x)
        cbind(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
    })
}
