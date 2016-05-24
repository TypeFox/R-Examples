cost_gaussian = function (x, data_x, data_y, prior_p, prior_st) 
{
    data_num = dim(data_x)[1]
    dm = dim(data_x)[2]
    hn = data_num^(-1/(4 + dm))
    tau2 = h = vector(, dm)
    for (k in 1:dm) {
        tau2[k] = exp(x[k])
        h[k] = sqrt(tau2[k]) * hn
    }
    hprod = prod(h)
    cont = exp(-0.5 * dm * log(2 * pi))
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
    logp = -0.5 * (data_num + prior_p) * log(0.5 * cv + 0.5 * prior_st)
    for (i in 1:dm) {
        logp = logp + x[i] + logpriorh2(tau2[i] * hn * hn) + log(hn * hn)
    }
    return(-logp)
}
