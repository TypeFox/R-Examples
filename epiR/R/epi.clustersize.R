"epi.clustersize" <- function(p, b, rho, epsilon.r, conf.level = 0.95){
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    D <- rho * (b - 1) + 1
    s <- p * epsilon.r
    numerator <- p * (1 - p) * D * z^2
    denominator <- (s^2 * b)
    clusters <- round((numerator/denominator) + 1, digits = 0)
    units <- clusters * b
    rval <- list(clusters = clusters, units = units, design = D)
    return(rval)
    }
