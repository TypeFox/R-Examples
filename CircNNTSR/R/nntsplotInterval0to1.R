nntsplotInterval0to1 <-
function (cpars = 1/sqrt(2*pi), M = 0, ...) 
{
    if (M == 0) {
        x <- rep(1, 2)
        return(plot(c(0, 1), x, type = "l", xlab = "theta"))
    }
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")

    nntsplotint <- function(theta) {
        res <- nntsDensityInterval0to1(theta, cpars, M)
	return(res)
    }
    return(curve(nntsplotint, 0, 1, xlab = "theta", ...))
}

