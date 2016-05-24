cost2_gaussian <- function(x, data_x, data_y, prior_st)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    tau2 = h = vector(,dm)
    hn = data_num^(-1/(4+dm))
    for(k in 1:dm)
    {
        tau2[k] = exp(x[k])
        h[k] = sqrt(tau2[k]) * hn
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
    cv = sum(cv_int) + prior_st
    return(0.5*cv)
}

