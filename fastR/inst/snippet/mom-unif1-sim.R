simulate <- function(size) {
    rdata <- runif(size)
    2 * mean(rdata) < max(rdata)
    }
mean(replicate(1000, simulate(6)))
mean(replicate(1000, simulate(12)))
mean(replicate(1000, simulate(24)))
