d2exp <- function (x, rate = 1, shift = 0, log = FALSE) 
{
    if (rate <= 0) {
        stop(paste("Rate must be larger than 0!", "\n"))
    }
    p <- exp(-(x - shift)/rate)/rate
    if (log) 
        p <- log(p)
    p
}


p2exp <- function (q, rate = 1, shift = 0, lower.tail = TRUE, log.p = FALSE) 
{
    if (rate <= 0) {
        stop(paste("Rate must be larger than 0!", "\n"))
    }
    p <- 1 - exp(-(q - shift)/rate)
    if (lower.tail == FALSE) 
        p <- 1 - p
    if (log.p) 
        p <- log(p)
    p
}


q2exp <- function (p, rate = 1, shift = 0, lower.tail = TRUE, log.p = FALSE) 
{
    if (rate <= 0) {
        stop(paste("Rate must be larger than 0!", "\n"))
    }
    if (log.p) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    x <- shift - rate * log(1 - p)
    x
}


r2exp <- function (n, rate = 1, shift = 0) 
{
    if (rate <= 0) {
        stop(paste("Rate must be larger than 0!", "\n"))
    }
    q2exp(runif(n), rate = rate, shift = shift)
}


