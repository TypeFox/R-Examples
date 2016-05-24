kern <- function(h, data_x, data_y, xm)
{
    data_num = nrow(data_x)
    dm = ncol(data_x)
    hn = data_num^(-1/(4 + dm))
    h = h[1:dm] * hn
    hprod = prod(h)
    cont = exp(-0.5 * dm * log(2.0 * pi))
    suma = sumb = 0
    for(i in 1:data_num)
    {
        suma = suma + data_y[i]
        sumb = sumb + data_y[i]^2
    }
    suma = suma/data_num
    sst = sumb - suma^2 * data_num
    sse = ase = 0
    for(i in 1:data_num)
    {
        suma = sumb = 0
        for(j in 1:data_num)
        {
            temp = ((data_x[i,] - data_x[j,])/h)^2
            weight = cont * exp(-0.5 * sum(temp))/hprod
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        sse = sse + (data_y[i] - suma/sumb)^2
        ase = ase + (xm[i] - suma/sumb)^2
    }
    r2 = 1 - sse/sst
    return(list(r2 = r2, mse = ase/data_num))
}

