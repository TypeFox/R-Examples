acvfFD <-
function(d, maxlag)
{
    x <- numeric(maxlag + 1)
    x[1] <- gamma(1 - 2 * d)/gamma(1 - d)^2
    for(i in 1:maxlag) 
        x[i + 1] <- ((i - 1 + d)/(i - d)) * x[i]
    x
}
