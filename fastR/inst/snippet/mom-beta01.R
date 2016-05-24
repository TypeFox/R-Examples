beta.mom <- function(x,lower=0.01,upper=100) {
    x.bar <- mean (x)
    n <- length(x)
    v <- var(x) * (n-1) / n
    R <- 1/x.bar - 1

    f <- function(a){             # note: undefined when a=0
        R * a^2 / ( (a/x.bar)^2 * (a/x.bar + 1) ) - v
    }

    u <- uniroot(f,c(lower,upper))

    return( c(shape1=u$root, shape2=u$root * R) )
}
x <- rbeta(50,2,5); beta.mom(x)
