mipvalue <- function (pvec)
{
    m <- length(pvec)
    pvec[pvec < .Machine$double.neg.eps] <- .Machine$double.neg.eps
    pvec[pvec==1] <- 1 - .Machine$double.neg.eps
    zi <- qnorm(1 - pvec)
    zm <- mean(zi)
    B <- var(zi)
    W <- 1
    T <- W + (1 + 1/m) * B
    nu <- (m - 1) * (1 + m * W/((m + 1) * B))^2
    mip <- pt(zm, df = nu, lower.tail = FALSE)
    return(mip)
}
