KurtosisAzimuthCoefficient <- function (azimuths) 
{
    n_elements = length(azimuths)
    sin_ = AddSinVectors(azimuths)
    cos_ = AddCosVectors(azimuths)
    CSub2 = (1/n_elements) * sum(cos(2 * ToRadians(azimuths)))
    SSub2 = (1/n_elements) * sum(sin(2 * ToRadians(azimuths)))
    RSub2 = ((CSub2^2 + SSub2^2)^0.5)
    R = MeanModule(azimuths)
    azimuthSub2 = atan(SSub2/CSub2)
    if ((SSub2 > 0) && (CSub2 > 0)) {
    }
    if ((SSub2 > 0) && (CSub2 < 0)) {
        azimuthSub2 = azimuthSub2 + pi
    }
    if ((SSub2 < 0) && (CSub2 > 0)) {
        azimuthSub2 = azimuthSub2 + 2 * pi
    }
    if ((SSub2 < 0) && (CSub2 < 0)) {
        azimuthSub2 = azimuthSub2 + pi
    }
    result = (RSub2 * (cos((azimuthSub2) - 2 * (atan(sin_/cos_)))) - 
        R^4)/((1 - R)^2)
    return(result)
}
