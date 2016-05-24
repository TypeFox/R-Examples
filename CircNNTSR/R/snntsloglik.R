snntsloglik <- function (data, cpars = 1, M = c(0,0)) 
{
    size <- length(cpars)
    if (size != (M[1] + 1)*(M[2] + 1)) 
        return("Length of cpars must be equal to (M[1]+1)*(M[2] + 1)")
    if (abs(sum(Mod(cpars)^2) - 1) > 1e-10) 
        return("sum of the squared norms of componentes greater than condition")
    auxcond1 <- sum(data[,1] > 2*pi) + sum(data[,1] < 0)
    auxcond2 <- sum(data[,2] > pi) + sum(data[,2] < 0)
    if (auxcond1>0)
	return("First column of the data matrix must have values between 0 and 2*pi")
    if (auxcond2>0)
	return("Second column of the data matrix must have values between 0 and pi")

    y <- sum(log(snntsdensity(data, cpars, M)))
    return(y)
}