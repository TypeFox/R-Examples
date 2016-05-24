CircularVariance <- function (azimuths) 
{
    module = MeanModule(azimuths)
    variance = 1.0 - module
    return(variance)
}
