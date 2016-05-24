von_mises_sim <-
function (n = 5000, k = 1, c = 0.3, noise = 0.2) 
{

    #convert c to radians
    c <- c * 2 * pi

    m <- optimise(von_mises_density, interval = c(0, 2 * pi), c = c, kappa = k, maximum = T)$objective
    x <- NULL


    noisenr <- length(which(runif(n, 0, 1) <= noise))
    signalnr <- n - noisenr
    
    while(length(x) < signalnr)
    {
          y <- runif(signalnr * m, c, 2 * pi - c)
          x <- c(x, y[runif(signalnr * m) * m < von_mises_density(y, c, k)])
    }

    x <- c(x[1 : signalnr], runif(noisenr, 0, 2 * pi))
    x <- sort(x) / (2 * pi)
    return(x)
}
