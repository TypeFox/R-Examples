nntsABDensity <-
function (theta, cpars = 1/sqrt(2*pi), M = 0) 
{
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norm of componentes greater than condition")

    ab <- nntsABcoefficients(cpars, M)
    y <- 1/(2 * pi)
    for (k in 1:M) {
        y <- y + ab[k] * cos(k * theta) + ab[(k + M)] * sin(k * 
            theta)
    }
    return(y)
}

