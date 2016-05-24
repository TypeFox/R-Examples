q.mle.norm.estimate <- function(complete.lifespans, censored.lifespans)
{
    mle(function(mean, sd, Q)
        {
            qsurvival.nllik("norm", complete.lifespans, censored.lifespans, Q, mean, sd)
        },
        method='L-BFGS-B',
        lower=list(mean=min(complete.lifespans),  sd=0.01, Q=0.000001),
        upper=list(mean=max(complete.lifespans),  
                   sd=(max(complete.lifespans) - min(complete.lifespans)),
                   Q=0.999999),
        start=list(mean=mean(complete.lifespans), 
                   sd=sd(complete.lifespans),
                   Q = 0.5))
}