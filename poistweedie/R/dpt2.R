dpt2 <-
function (p, n, mu, lambda, theta0) 
{
    y <- double((n + 1))
    dpt <- double((n + 1))
    w <- double(1)
    K <- double(1)
    c0 <- double(1)
    suma <- double(n)
    sumv <- matrix(0, n, n)
	## pr <- double(n)
   
    ## fixons la précision a 14
	## pr<- 14
    w <- omega(p,mu,theta0)
    y <- c(0:n)
    for (i in 1:(n + 1)) {
        if (theta0 == 0) {
            dpt[i] <- gam1.1(i - 1, lambda) * exp(w * (i - 1) + 
                lambda * log(1 - exp(w)))
        }
        else {
            dpt[i] <- gam1.1(i - 1, lambda) * exp((i - 1) * (w + 
                log(1/(1 - theta0))) + lambda * log((exp(w) - 
                1 + theta0)/(theta0 - 1)))
        }
    }
    dpt
}

