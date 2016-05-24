MeanAngle <- function (angles) 
{
    sin_ = AddSinVectors3D(angles)
    cos_ = AddCosVectors3D(angles)
    angle = atan(sin_/cos_)
    angle = ToSexagesimal3D(angle)
    return(angle)
}
