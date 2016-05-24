logdensity_admkr = function(tau2, cpost)
{
    dm = ncol(cpost) - 1
    len = nrow(cpost)
    band = vector(, dm + 1)
    for(j in 1:(dm+1))
    {
        temp = tem2 = 0
        for(i in 1:len)
        {
            temp = temp + cpost[i,j]
            tem2 = tem2 + cpost[i,j]^2
        }
        sigma = sqrt(tem2/len - (temp/len)^2)
        temp = exp(1.0/(dm + 5) * log(4/(dm + 3)))
        band[j] = temp * sigma * exp(-1/(dm + 5) * log(len))
    }
    hprod = prod(band)
    cont = exp(-0.5 * (dm+1) * log(2.0 * pi))
    xsum = 0
    for(i in 1:len)
    {
        temp = 0
        for(j in 1:(dm + 1))
        {
            tem2 = (tau2[j] - cpost[i,j])/band[j]
            temp = temp + tem2^2
        }
        xsum = xsum + cont * exp(-0.5 * temp)/hprod
    }
    hatf = log(xsum/len)
    return(hatf)
}

logdensity_admkr = function(tau2, cpost)
{
    dm = ncol(cpost)
    len = nrow(cpost)
    band = vector(,dm)
    for(j in 1:dm)
    {
        temp = tem2 = 0
        for(i in 1:len)
        {
            temp = temp + cpost[i,j]
            tem2 = tem2 + cpost[i,j]^2
        }
        sigma = sqrt(tem2/len - (temp/len)^2)
        temp = exp(1.0/(dm+5) * log(4/(dm+3)))
        band[j] = temp * sigma * exp(-1/(dm+5) * log(len))
    }
    hprod = prod(band)
    cont = exp(-0.5 * (dm+1))
}
