q.mle.emg.estimate <- function(complete.lifespans, censored.lifespans)
{
  mle(function(mu, sigma, lambda, Q)
        {
           qsurvival.nllik("emg", complete.lifespans, censored.lifespans, Q, mu, sigma, lambda)
        },
        method='L-BFGS-B',
        lower=list(mu=8,  sigma=0.1, lambda=0.01, Q=0),
        upper=list(mu=30, sigma=(max(complete.lifespans) - min(complete.lifespans)),
                   lambda=10,    Q=0.999999),
        start=list(mu=median(complete.lifespans), 
                   sigma=sd(complete.lifespans),
                   lambda=1/mean(complete.lifespans),
                   Q = 0.5))
}