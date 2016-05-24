gsHP <-
function (ystar, rstar, distr) 
{
    if (distr == 1) {
        mu0 <- 0
        s0 <- 0.01
        q1 <- 0.1
        q2 <- 0.1
        a <- q1 + rstar/2
        if (rstar > 1) {
            b <- q2 + (rstar - 1) * var(ystar)/2 + s0 * rstar * 
                (mean(ystar) - mu0)^2/2/(s0 + rstar)
        }
        else {
            b <- q2 + s0 * rstar * (mean(ystar) - mu0)^2/2/(s0 + 
                rstar)
        }
        t2 <- rgamma(1, shape = a, rate = b)
        a <- (s0 * mu0 + sum(ystar))/(s0 + rstar)
        b <- (s0 + rstar) * t2
        t1 <- rnorm(1, mean = a, sd = 1/sqrt(b))
        mu.py0 <- t1
        sigma.py0 <- 1/sqrt(t2)
    }
    else if (distr == 2) {
        q1 <- 0.01
        q2 <- 0.01
        t1 <- rgamma(1, shape = q1 + rstar, rate = q2 + sum(ystar))
        mu.py0 <- sigma.py0 <- 1/t1
    }
    else if (distr == 3) {
        q1 <- 0.01
        q2 <- 0.01
        t1 <- rgamma(1, shape = q1 + rstar, rate = q2 - sum(log(ystar)))
        mu.py0 <- t1/(t1 + 1)
        sigma.py0 <- sqrt(t1/(t1 + 1)^2/(t1 + 2))
    }
    else {
        stop("Argument \"distr\" should be defined numeric with possible values 1,2 or 3")
    }
    return(list(mu.py0 = mu.py0, sigma.py0 = sigma.py0))
}
