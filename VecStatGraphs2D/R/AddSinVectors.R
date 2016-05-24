AddSinVectors <- function (vectors) 
{
    sin_sum = 0
    radians = ToRadians(vectors)
    radians_sin = sin(radians)
    sin_sum = sum(radians_sin)
    return(sin_sum)
}
