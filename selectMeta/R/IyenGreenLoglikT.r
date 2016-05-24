IyenGreenLoglikT <- function(para, t, q, N, type = 1){

    theta <- para[1]
    b <- para[2]

    iG <- 0
    int <- 0
    
    for (i in 1:length(q)){
        iG <- iG + log(IyenGreenWeight(t[i], b, q[i], type = type, alpha = 0.05))
        tmp <- integrate(f = normalizeT, lower = -Inf, upper = Inf, stop.on.error = FALSE, theta = theta, b = b, q = q[i], N = N[i], type = type)$value
        int <- int + log(tmp)
        }

    loglik <- suppressWarnings(sum(dt(t, ncp = sqrt(N / 2) * theta, df = q, log = TRUE)) + iG - int)
    return(loglik)
}
