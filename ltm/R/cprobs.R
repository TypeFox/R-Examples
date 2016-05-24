cprobs <-
function (betas, z, eps = .Machine$double.eps^(1/3)) {
    lapply(betas, function (x, z) {
        nx <- length(x)
        out <- plogis(x[-nx] - matrix(x[nx] * z, nx - 1, length(z), TRUE))
        if (any(ind <- out == 1))
            out[ind] <- 1 - eps
        if (any(ind <- out == 0))
            out[ind] <- eps
        rbind(out, 1)        
    }, z = z)
}
