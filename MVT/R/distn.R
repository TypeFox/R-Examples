dweights <-
function(x, eta = .25, dim, log = FALSE, scaled = TRUE)
{
    log.density <- log
    if (scaled) {
        y <- dbeta(x, .5 / eta, dim / 2, 0, log)
    }
    else {
        max <- (1 + eta * dim) / (1 - 2 * eta)
        a <- .5 / eta
        b <- dim / 2
        y <- -a * (eta * (dim - 2) + 1) * log(max) + (a - 1) * log(x) + (b - 1) * log(max - x) - lbeta(a, b)
        y <- ifelse(log.density, y, exp(y))
    }
    y
}

pweights <-
function(q, eta = .25, dim, lower.tail = TRUE, log.p = FALSE, scaled = TRUE)
{
    if (scaled) {
        y <- pbeta(q, .5 / eta, dim / 2, 0, lower.tail, log.p)
    }
    else {
        max <- (1 + eta * dim) / (1 - 2 * eta)
        y <- pbeta(q / max, 5. / eta, dim / 2, 0, lower.tail, log.p)
    }
    y
}

qweights <-
function(p, eta = .25, dim, lower.tail = TRUE, log.p = FALSE, scaled = TRUE)
{
    y <- qbeta(p, .5 / eta, dim / 2, 0, lower.tail, log.p)
    if (!scaled) {
        max <- (1 + eta * dim) / (1 - 2 * eta)
        y <- max * y
    }
    y
}
