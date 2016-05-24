nicerange = function (x, beta = 0.1)
{
    ab <- range(x)
    del <- ((ab[2] - ab[1]) * beta)/2
    return(c(ab + c(-del, del)))
}
