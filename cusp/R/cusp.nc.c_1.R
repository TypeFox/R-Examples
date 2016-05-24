`cusp.nc.c` <-
function (alpha, beta, ..., keep.order = TRUE) 
{
    i = order(alpha, beta)
    nc <- cusp.nc.C(alpha[i], beta[i], ...)
    if (keep.order) 
        nc$value[(1:length(alpha))[i]]
    else nc$value
}

