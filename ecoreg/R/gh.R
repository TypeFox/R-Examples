
## gauss.hermite() and hermite() functions by Jim Lindsey from his rmutil package

gauss.hermite <- function (points, iterlim = 50) 
{
    x <- w <- rep(0, points)
    m <- (points + 1)/2
    for (i in 1:m) {
        z <- if (i == 1) 
            sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
        else if (i == 2) 
            z - sqrt(points)/z
        else if (i == 3 || i == 4) 
            1.9 * z - 0.9 * x[i - 2]
        else 2 * z - x[i - 2]
        for (j in 1:iterlim) {
            z1 <- z
            p <- hermite(points, z)
            z <- z1 - p[1]/p[2]
            if (abs(z - z1) <= 1e-15) 
                break
        }
        if (j == iterlim) 
            warning("iteration limit exceeded")
        x[points + 1 - i] <- -(x[i] <- z)
        w[i] <- w[points + 1 - i] <- 2/p[2]^2
    }
    r <- cbind(x * sqrt(2), w/sum(w))
    colnames(r) <- c("Points", "Weights")
    r
}

hermite <- function (points, z) 
{
    p1 <- 1/pi^0.4
    p2 <- 0
    for (j in 1:points) {
        p3 <- p2
        p2 <- p1
        p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
    }
    pp <- sqrt(2 * points) * p2
    c(p1, pp)
}

## Use Gauss-Hermite approximation to integrate the vectorised function h over its first vector argument. 
## CJ. Uses gauss.hermite() function from Jim Lindsey's rmutil package

integrate.gh <- function(h, n=1, points=10, mu=0, scale=1, ...){
    gh <- gauss.hermite(points)
    pts <- gh[,"Points"]
    wts <- gh[,"Weights"]
    s <- 0
    h2 <- function(x, ...) h(x, ...) / dnorm(x, mu, scale)
    for (i in 1:points) {
        s <- s + (wts[i])*h2(mu + scale*pts[i], ...) 
    }
    s
}
