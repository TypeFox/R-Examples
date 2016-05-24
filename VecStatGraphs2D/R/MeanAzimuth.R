MeanAzimuth <- function (azimuths) 
{
    sin_ = AddSinVectors(azimuths)
    cos_ = AddCosVectors(azimuths)
    azimuth = atan(sin_/cos_)
    azimuth = ToSexagesimal(azimuth)
    if ((sin_ > 0) && (cos_ > 0)) {
        azimuth = azimuth
    }
    if ((sin_ > 0) && (cos_ < 0)) {
        azimuth = azimuth + 180
    }
    if ((sin_ < 0) && (cos_ > 0)) {
        azimuth = azimuth + 360
    }
    if ((sin_ < 0) && (cos_ < 0)) {
        azimuth = azimuth + 180
    }
    return(azimuth)
}
