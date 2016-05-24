logpriors_admkr <- function(h2, prior_alpha, prior_beta)
{
    logf = 0
    dm = length(h2)
    for(i in 1:dm)
    {
        logf = logf + logpriorh2(h2[i], prior_alpha = prior_alpha, prior_beta = prior_beta)
    }
    return(logf)
}
