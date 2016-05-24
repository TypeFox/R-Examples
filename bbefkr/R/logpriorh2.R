logpriorh2 <- function(h2, prior_alpha, prior_beta)
{
    logp = prior_alpha * log(prior_beta) - lgamma(prior_alpha) - 1.0*(prior_alpha + 1.0)*log(h2) - prior_beta/h2	
    return(logp)
}
