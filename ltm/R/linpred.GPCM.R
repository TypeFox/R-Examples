linpred.GPCM <-
function (betas, z, IRT.param = TRUE) {
    lapply(betas, function (x) {
        nx <- length(x)
        if (IRT.param)
            t(x[nx] * outer(z, x[-nx], "-"))
        else
            outer(x[-nx], x[nx] * z , "+")
    })
}
