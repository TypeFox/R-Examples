CircularStandardDeviation <- function (azimuths) 
{
    variance = CircularVariance(azimuths)
    deviation = (-2 * log(1 - variance))^0.5
    return(deviation)
}
