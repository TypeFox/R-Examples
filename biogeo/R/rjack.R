rjack <-
function (d) 
{
    xx <- d
    d <- unique(d)
    rng <- diff(range(d))
    mx <- mean(d)
    n <- length(d)
    n1 <- n - 1
    t1 <- (0.95 * sqrt(n)) + 0.2
    x <- sort(d)
    y <- rep(0, n1)
    for (i in 1:n1) {
        x1 <- x[i + 1]
        if (x[i] < mx) {
            y[i] <- (x1 - x[i]) * (mx - x[i])
        }
        else {
            y[i] <- (x1 - x[i]) * (x1 - mx)
        }
    }
    my <- mean(y)
    z <- y/(sqrt(sum((y - my)^2)/n1))
    out <- rep(0, length(xx))
    if (any(z > t1)) {
        f <- which(z > t1)
        v <- x[f]
        if (v < median(x)) {
            xa <- (xx <= v) * 1
            out <- out + xa
        }
        if (v > median(x)) {
            xb <- (xx >= v) * 1
            out <- out + xb
        }
    }
    else {
        out <- out
    }
    f <- which(out == 1)
}
