AddSinVectors3D <- function (vectors) 
{
    sin_sum = 0
    h = 1
    radians = ToRadians3D(vectors)
    radians_sin = sin(radians)
    sin_sum = sum(radians_sin)
    return(sin_sum)
}
