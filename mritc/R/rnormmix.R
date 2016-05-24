rnormmix <- function(n, prop, mu, sigma){

    k <- length(prop)
    if (!(k == length(mu) && k == length(sigma)))
        stop("The dimensions of 'prop', 'mu' and 'sigma' do NOT match.")
    if (!(all(prop > 0) && all(prop < 1)))
        stop("'prop' has to be between 0 and 1")
    if (! isTRUE(all.equal(sum(prop), 1)))
        stop("Sum of 'prop' has to be 1")
    if (!all(sigma > 0))
        stop("All 'sigma's have to be positive")

    c <- sample(1:k, size=n, replace=TRUE, prob=prop)
    y <- rep(0, n)
    for(i in 1:k){
        ci <- c==i
        y[ci] <- rnorm(sum(ci), mean=mu[i], sd=sigma[i])
    }
    cbind(y,c)
}
