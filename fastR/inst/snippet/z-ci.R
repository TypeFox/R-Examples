ci <- function(x, sd=100, conf.level=0.95) { 
    alpha = 1 - conf.level
    n = length(x)
    zstar <- - qnorm(alpha/2)
    interval <- mean(x)  + c(-1,1) * zstar * sd / sqrt(n)
    return(list(conf.int=interval, estimate=mean(x)))
    }
