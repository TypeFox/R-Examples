DearBeggToMinimize <- function(vec, y, u, lam){
    n <- length(y)
    k <- 1 + floor(n / 2)
    w <- vec[1:k]
    teststat <- abs(y) / u

    res <- Inf
    crit1 <- min(diff(w)) >= 0
    crit2 <- max(w) > 0

    if (crit1 & crit2){
        theta <- vec[k + 1]
        sigma <- vec[k + 2]
        hij <- Hij(theta = theta, sigma = sigma, y = y, u = u, teststat = teststat)$Hij
        res <- -DearBeggLoglik(w = w, theta = theta, sigma = sigma, y = y, u = u, hij = hij, lam = lam)$LL
        }
    
    return(res)
    }
