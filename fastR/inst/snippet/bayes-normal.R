x <- c(20,24,27,28,28,28,29,30,30,30,30,32,33,34,35,38)
mean(x)
sd(x)
posterior <- function(x,mu0,sigma0,sigma=5) {
    n <- length(x)
    N <- (n*mean(x)/sigma^2 + mu0/sigma0^2)
    D <- (n/sigma^2 + 1/sigma0^2)
    mu1 <- N/D; sigma1 <- sqrt(1/D)  
    precision1 <- D
    precision0 <- 1/sigma0^2
    precision.data <- n/sigma^2
    return(cbind(mu1,sigma1,precision1,precision0,precision.data))
    }
posterior(x,20,1)
posterior(x,20,4)
posterior(x,20,16)
posterior(x,20,1000)
5/sqrt(length(x))
