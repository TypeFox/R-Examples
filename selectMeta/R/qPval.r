qPval <- function(prob, u, theta, sigma2){

    tmp.fun <- function(q, prob, u, theta, sigma2){return(pPval(q, u, theta, sigma2) - prob)}

    if (tmp.fun(q = 10^-10, prob, u, theta, sigma2) >= 0){res <- 10^-5} else {res <- uniroot(tmp.fun, interval = c(10^-10, 1), prob, u, theta, sigma2, tol = 10^-8)$root}
    return(res)
}

