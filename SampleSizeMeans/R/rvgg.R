`rvgg` <-
function(alpha, beta, n)
{
# Generation of Gamma-gamma variables with parameters given in the vectors
#   alpha, beta and n
#
# If any of the parameters doesn't have the same length as the longest,
# it is replicated a bunch of times to be the same length as others
        m1 <- length(alpha)
        m2 <- length(beta)
        m3 <- length(n)
        m <- max(m1, m2, m3)
        alpha <- rep(alpha, m/m1)
        beta <- rep(beta, m/m2)
        n <- rep(n, m/m3)
        l <- 1/beta * rgamma(m, alpha)
        1/l * rgamma(m, n)
}

