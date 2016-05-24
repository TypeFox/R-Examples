dpt1 <-
function(p, n, mu, lambda, theta0) 
{
    k <- numeric(1)
    i <- numeric(1)
    a <- matrix(0, n, n)
    fac <- double(1)
    w <- double(1)
	## pr <- double(1)
    alpha <- double(1)
    K <- double(1)
    suma <- matrix(0, n, n)
    dpt <- double((n + 1))
    
    ## fixons la pr\'{e}cision a 14
	## pr<- 14
    w <- omega(p,mu,theta0)
    K <- exp(exp(w) - 1 + theta0) - exp(theta0)
    c0 <- exp(lambda * (exp(theta0 - 1)) * (1 - exp(1)))
    if (n >= 1) {
        a[1, 1] <- exp(w + (theta0 - 1 + log(lambda)) - lambda * 
            K)
        if (n >= 2) {
            a[2, 1] <- a[1, 1] * exp(w - log(2))
            a[2, 2] <- a[1, 1] * exp(w - log(2) + (theta0 - 1 + 
                log(lambda)))
            if (n > 2) {
                for (i in 3:n) {
                  a[i, 1] <- a[(i - 1), 1] * exp(w - log(i))
                  a[i, i] <- a[i - 1, i - 1] * exp(w - log(i) + 
                    (theta0 - 1 + log(lambda)))
                  for (k in 2:(i - 1)) {
                    a[i, k] <- exp(w - log(i) + log(exp(theta0 - 
                      1 + log(lambda)) * a[i - 1, k - 1] + k * 
                      a[i - 1, k]))
                  }
                }
            }
            for (i in 1:n) {
                suma[i] <- sum(a[i, ])
            }
        }
    }
    dpt[1] <- c0 * exp(-lambda * K)
    if (n >= 1) {
        for (i in 2:(n + 1)) {
            dpt[i] <- c0 * suma[i - 1]
        }
    }
    dpt
}

