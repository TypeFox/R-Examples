CircularDispersal <- function (azimuths) 
{
    n = length(azimuths)
    CSub2 = (1/n) * sum(cos(2 * ToRadians(azimuths)))
    SSub2 = (1/n) * sum(sin(2 * ToRadians(azimuths)))
    RSub2 = ((CSub2^2 + SSub2^2)^0.5)
    R = MeanModule(azimuths)
    result = (1 - RSub2)/(2 * (R^2))
    return(result)
}
