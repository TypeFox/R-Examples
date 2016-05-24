MeanModule <- function (azimuths) 
{
    n = length(azimuths)
    sin_ = AddSinVectors(azimuths)
    cos_ = AddCosVectors(azimuths)
    module = sqrt((sin_ * sin_) + (cos_ * cos_))
    mean_module = module / n
    return(mean_module)
}
