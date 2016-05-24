euclidian <-
function (v1, v2)
{
    dv <- v1 - v2
    return(sqrt(sum(dv^2)))
}
