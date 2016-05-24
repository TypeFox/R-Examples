dptp <-
function(p, n, mu, lambda, theta0) 
{
    k <- numeric(1)
    para <- double(1)
    i <- numeric(1)
    a <- matrix(0, n, n)
    sumv <- matrix(0, n, n)
    KCsum <- double(n)
    y <- double((n + 1))
    dpt <- double((n + 1))
    w <- double(1)
    K1 <- double(1)
    c0 <- double(1)
    alpha <- double(1)
	## pr <- double(n)
   
    ## fixons la précision a 14
	##pr<- 14
    
    alpha <- (p - 2)/(p - 1)
    cons <- log((1 - theta0)/(1 - alpha))
    para <- (log(lambda) + alpha * log((1 - theta0)/(1 - alpha)))
    w <- omega(p,mu,theta0)
    if (theta0 == 0) {
        K1 <- (1/(2 - p)) * (((1 - p) * (exp(w) - 1))^(alpha))
    }
    else {
        K1 <- (1/(2 - p)) * (((1 - p) * (exp(w) - 1 + theta0))^(alpha) - 
            ((1 - p) * theta0)^(alpha))
    }
    if (theta0 == 0) {
        c0 <- exp(((lambda * (alpha - 1))/alpha) * ((1/(1 - alpha))^(alpha)))
    }
    else {
        c0 <- exp(((lambda * (alpha - 1))/alpha) * ((((1 - theta0)/(1 - 
            alpha))^(alpha)) - (((-theta0)/(1 - alpha))^(alpha))))
    }
    if (n >= 1) {
        a <- diag(1, n)
        a[1, 1] <- exp(w - cons + para - lambda * K1)
        if (n >= 2) {
            a[2, 1] <- exp(2 * (w - log(1 - theta0)) + 2 * log(1 - 
                alpha) + para - lambda * K1)/2
            a[2, 2] <- a[1, 1] * exp(w - cons - log(2) + para)
            if (n >= 3) {
                for (i in 3:n) {
                  a[i, 1] <- exp(i * (w - log(1 - theta0)) + 
                    2 * log(1 - alpha) + para - lambda * K1 + 
                    log((gam1.1(i, -alpha)/((2 * gam1.1(2, -alpha))))))
                  a[i, i] <- a[i - 1, i - 1] * exp(w - cons - 
                    log(i) + para)
                  for (k in 2:(i - 1)) {
                    a[i, k] <- exp((w - log(i) - cons) + log((a[i - 
                      1, k - 1] * exp(para) + ((i - 1 - k * alpha)/(1 - 
                      alpha)) * a[i - 1, k])))
                  }
                }
            }
        }
        for (i in 1:n) {
            KCsum[i] <- sum(a[i, ])
        }
    }
    dpt[1] <- c0 * exp(-lambda * K1)
    if (n >= 1) {
        for (i in 2:(n + 1)) {
            dpt[i] <- c0 * KCsum[i - 1]
        }
    }
    dpt
}

