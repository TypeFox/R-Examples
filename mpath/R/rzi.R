### random number generation for zero-inflated response variable
### x is the design matrix of count model
### z is the design matrix of zero model
### a is the coefficient for x
### b is the coefficient for z
rzi <- function(n, x, z, a, b, theta=1, family=c("poisson", "negbin", "geometric"), infl=TRUE){
    family <- match.arg(family)
    if(family=="geometric" && theta!=1) stop("theta should be 1 for family='geometric'")
    mu <- exp(a[1] + x %*% a[-1])
    p <- exp(b[1] + z %*% b[-1])
    w <- 1/(1+1/p)  ### logit link
    if(family=="poisson")
        y <- rpois(n, lambda=mu)
    else y <- rnegbin(n, mu=mu, theta=theta) ###check
    if(infl){
        y[runif(n) < w] <- 0
        cat("zero inflation", mean(w), "\n")
        cat("proportion of y=0", (length(y[y==0])/length(y)), "\n")
    }
    return(y)
}

