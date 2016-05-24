infoprobs <-
function (betas, z) {
    cpr <- cprobs(betas, z)
    ipr <- iprobs(betas, z)
    sum.cprs <- lapply(cpr, function (x) {
        nr <- nrow(x)
        t((1 - rbind(x[1, ], x[-nr, ] + x[-1, ]))^2)
    })
    betas. <- sapply(betas, function (x) x[length(x)])
    for (i in 1:length(betas))
        sum.cprs[[i]] <- betas.[i]^2 * ipr[[i]] * sum.cprs[[i]]
    do.call(cbind, lapply(sum.cprs, rowSums))
}
