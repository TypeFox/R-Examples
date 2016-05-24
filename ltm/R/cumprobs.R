cumprobs <-
function (betas, z,  lower = TRUE) {
    lapply(betas, function (x, z) {
        nx <- length(x)
        out <- plogis(x[-nx] - matrix(x[nx] * z, nx - 1, length(z), TRUE))
        if (lower) t(rbind(out)) else 1 - t(rbind(out))
    }, z = z)
}
