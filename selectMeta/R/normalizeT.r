normalizeT <- function(s, theta, b, q, N, type = 1, alpha = 0.05){

    p1 <- suppressWarnings(dt(s, ncp = sqrt(N / 2) * theta, df = q))
    p2 <- IyenGreenWeight(s, b, q, type, alpha)
    res <- p1 * p2
    return(res)
    }
