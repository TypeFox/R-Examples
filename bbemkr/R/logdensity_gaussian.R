logdensity_gaussian <- function(tau2, cpost)
{
    dm = ncol(cpost) - 1
    length = nrow(cpost)
    band = vector(, dm + 1)
    for(j in 1:(dm + 1))
    {
        temp = tem2 = 0
        for(i in 1:length)
        {
            temp = temp + cpost[i,j]
            tem2 = tem2 + cpost[i,j]^2
        }
        sigma = sqrt(tem2/length - (temp/length)^2)
        temp = exp(1/(dm + 5) * log(4/(dm + 3)))
        band[j] = temp * sigma * exp(-1/(dm + 5) * log(length))
    }
    hprod = prod(band)
    cont = exp(-0.5 * (dm + 1) * log(2 * pi))
    sum = 0
    for(i in 1:length)
    {
        temp = 0
        for(j in 1:(dm + 1))
        {
            temp = temp + ((tau2[j] - cpost[i,j])/band[j])^2
        }
        sum = sum + cont * exp(-0.5 * temp)/hprod
    }
    hatf = log(sum/length)
    return(hatf)
}

