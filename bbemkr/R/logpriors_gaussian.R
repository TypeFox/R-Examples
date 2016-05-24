logpriors_gaussian <- function(h2, data_x, prior_p, prior_st)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    hn = data_num^(-1/(4 + dm))
    sigma2 = h2[dm+1]
    logp = 0
    for(i in 1:dm)
    {
        logp = logp + logpriorh2(h2[i] * hn * hn) + log(hn * hn)
    }
    logp = logp + 0.5 * prior_p * log(0.5 * prior_st) - lgamma(0.5 * prior_p) - (0.5 * prior_p + 1.0) * log(sigma2) - 0.5 * prior_st/sigma2
    return(logp)
}

