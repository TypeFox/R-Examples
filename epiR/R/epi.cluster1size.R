"epi.cluster1size" <- function(n, mean, var, epsilon.r, method = "mean", conf.level = 0.95){
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    if (method == "total") {
        Vsq <- var / mean^2
        numerator <- z^2 * n * Vsq
        denominator <- z^2 * Vsq + (n - 1) * epsilon.r^2
        rval <- round(numerator/denominator, digits = 0)
        }
    
    if (method == "mean") {
        Vsq <- var / mean^2
        numerator <- z^2 * n * Vsq
        denominator <- z^2 * Vsq + (n - 1) * epsilon.r^2
        rval <- round(numerator/denominator, digits = 0)
        }
    
    if (method == "mean.per.unit") {
        Vsq <- var / mean^2
        numerator <- z^2 * n * Vsq
        denominator <- z^2 * Vsq + (n - 1) * epsilon.r^2
        rval <- round(numerator/denominator, digits = 0) 
       }
    if (method == "proportion") {
        if (length(var) != 2) 
           stop("Error: var must be of length 2")
        if (length(mean) != 2)
            stop("Error: mean must be of length 2")
        rval <- 'Not implemented yet!' 
       }
    return(rval)
}
