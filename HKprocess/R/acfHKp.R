acfHKp <- function(H,maxlag) {
    h2 <- 2*H
    k <- 1:maxlag
    c(1,0.5 * ((k + 1)^h2 - 2 * k^h2 + (k - 1)^h2))
}