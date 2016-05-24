estimRawResid <- function(U, Delta, covTerm)
## Estimates raw residuals (inverse normal scores) of a right-censored sample.
##
## 'U' is a vector of survival times, 'Delta' a censoring indicator, and
## 'covTerm' a scalar product of covariates with regression parameter estimates.
{
    tObs <- sort(U[Delta == 1])
    t1st <- tObs[1]
    t1stInd <- (U >= t1st)
    
    ## Estimating the step function of the transformation model with censored
    ## data.
    Hval <- numeric(length(tObs))
    Hval[1] <- uniroot(equationHtk, lower=-100, upper=100, tiInd=t1stInd,
                       covTerm=covTerm)$root

	for (i in 2:length(tObs))
	{
        ti <- tObs[i]
        tiInd <- (U >= ti)
        rskVec <- exp(covTerm + Hval[i - 1])
        Hval[i] <- Hval[i - 1] + 1 / sum(rskVec[tiInd])
	}

    Hfct <- stepfun(tObs, c(-Inf, Hval))
    rawResid <- qnorm(1 - exp(- exp(Hfct(U) + covTerm)))
    return(rawResid)
}