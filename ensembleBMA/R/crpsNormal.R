crpsNormal <-
function(sd, weights, biasCoefs, ensembleData)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  #a couple of helper functions
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  absExp <- function(mu, sig) 
  {
   (sqrt(2)* sig)*exp(-(mu/sig)^2/2)/sqrt(pi) + 
                       mu * erf((sqrt(2)*mu)/(2*sig))
  }


  # Expression of the CRPS formula and the E|x| if x ~ N(mu,sigma^2)

  # CRPS = .5 sum( sum( w(i)w(j) a( u(i) - u(j), sigma(i)^2 + sigma(j)^2) ) ) 
  #   - sum( w(i) a( mu(i) - obs, sigma(i)^2 )
  # here, a(u, sigma^2) is from E|X| with X ~ N(u, sigma^2)
  # Using Maple, I get Expected value of abs( X ) with X ~ N > >
  # (sigma*sqrt(2)*exp(-1/2/sigma^2*mu^2)+mu*erf(1/2/sigma*mu*2^(1/2))
  # *sqrt(Pi)) > / Pi^(1/2) > > 
  # where erf is the error function.
 
    nForecasts <- length(weights)

    if (length(sd) == 1) sd <- rep(sd, nForecasts)
    VAR <- sd*sd

    MEAN <- sweep(ensembleForecasts(ensembleData), MARGIN = 2, FUN = "*", 
                  STATS = biasCoefs[2,])
    MEAN <- sweep(MEAN, MARGIN = 2, FUN = "+", STATS = biasCoefs[1,])

    obs <- dataVerifObs(ensembleData)
    nObs <- length(obs)

    crpsTP <- numeric(nObs)

    for (l in 1:nObs) {

       M <- is.na(MEAN[l,])

       W <- weights/sum(weights[!M])

       crps1 <- crps2 <- 0

  # Begin computing the first term in the CRPS formula.  
  # This is a double sum since it is over w(i)*w(j) for all i and j.

       for (i in (1:nForecasts)[!M]) 
         {
            for (j in (1:nForecasts)[!M]) 
               {
                  tvar <- VAR[i] + VAR[j]  # total variance
                  tsd <- sqrt(tvar)          # total standard deviation
                  tmean <- MEAN[l,i] - MEAN[l,j]
                  temp <- absExp(tmean,tsd)
                  term <- (W[i]*W[j])*temp
                  crps1 <- crps1 + term
               }

             tvar <- VAR[i]              # total variance
             tsd <- sqrt(tvar)            # total standard deviation
             tmean <- MEAN[l,i] - obs[l]
             crps2 <- crps2 + W[i]*absExp(tmean,tsd)
        }

    # Using Szekely's expression for the CRPS, 
    # the first sum and second are put together to compute the CRPS.

     crpsTP[l]  <- crps2 - crps1/2     
    }

  mean(crpsTP)
}

