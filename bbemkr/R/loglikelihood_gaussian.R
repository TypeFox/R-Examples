loglikelihood_gaussian <- function(h2, data_x, data_y)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    hn = data_num^(-1/(4 + dm))
    h = vector(,dm)
    sigma2 = h2[dm + 1]
    for(i in 1:dm)
    {
        h[i] = sqrt(h2[i]) * hn
    }
    hprod = prod(h)
    cont = exp(-0.5 * dm * log(2.0 * pi))
    
    cv_int = vector(,data_num)
    for(i in 1:data_num)
    {
        temp = (sweep(data_x[-i,], 2, data_x[i,])/h)^2
        weight = cont * exp(-0.5 * apply(temp,1,sum))/hprod
        suma = sum(weight * data_y[-i])
        sumb = sum(weight)
        cv_int[i] = (data_y[i] - suma/sumb)^2
    }
    cv = sum(cv_int)
    logp = -0.5 * data_num * log(2.0 * pi * sigma2) - cv/(2.0 * sigma2)
    return(logp)
}

