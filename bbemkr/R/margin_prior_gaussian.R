margin_prior_gaussian <- function(band2, data_x, prior_p, prior_st)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    hn = data_num^(-1/(4 + dm))
    sigma2 = band2[1]
    logp = 0.5 * prior_p * log(0.5 * prior_st) - lgamma(0.5 * prior_p) -
            (0.5 * prior_p + 1) * log(sigma2) - 0.5 * prior_st/sigma2
    for(i in 2:(dm + 1))
    {
        logp = logp + logpriorh2(band2[i] * hn * hn) + log(hn * hn)
    }
    return(logp)
}

