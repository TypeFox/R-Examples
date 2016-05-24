nntsplot <-function (cpars = 1/sqrt(2*pi), M = 0, ...) 
{
    if (M == 0) {
        x <- rep(1/(2 * pi), 2)
        return(plot(c(0, 2 * pi), x, type = "l", ...))
     }
     size <- length(cpars)
     if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")

    nntsplotint <- function(theta) {
        res <- nntsdensity(theta, cpars, M)
        return(res)
    }
    return(curve(nntsplotint, 0, 2 * pi, ...))
}

