rPval <- function(n, u, theta, sigma2, seed = 1){

    ## generate a random number from the p-value distribution
    ## note: u is the trial-specific standard error, not a uniform random number!
    ## note that the sign of the effect sizes depends on theta!
    set.seed(seed)
    rn <- rep(NA, n)
    yn <- rep(NA, n)
    runi <- runif(n)

    for (i in 1:n){
        rn[i] <- qPval(runi[i], u, theta, sigma2)

        y <- - u * qnorm(rn[i] / 2)
        yn[i] <- y
        sign_i <- sign(runif(1, -1, 1))
        if (sign_i == -1){yn[i] <- 2 * theta - y}

        if (round(i / 100) == i / 100){print(paste("random number: ", i, " of ", n, " generated", sep = ""))}
        }    

    res <- list("rn" = rn, "yn" = yn)
    return(res)
}





