AddCosVectors <- function (vectors) 
{
    cos_sum = 0
    radians = ToRadians(vectors)
    radians_cos = cos(radians)
    cos_sum = sum(radians_cos)
    return(cos_sum)
}
