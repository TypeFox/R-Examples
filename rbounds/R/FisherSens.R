FisherSens <- function(totalN,
                       treatedN,
                       totalSuccesses,
                       treatedSuccesses,
                       Gammas) { 
  ## 'FisherSens' 
  ## Created by Devin Caughey on 1 March 2010
  ## Last modified: 4 December 2010
  ##
  ## This function performs a sensitivity analysis for Fisher's Exact Test
  ## for count (binary) data. It is derived from Section 4.4 of Paul
  ## Rosenbaum’s ``Observational Studies" (2nd Ed., 2002). The test is
  ## one-sided; that is, the alternative hypothesis to the null is
  ## that the number of "successes" (however defined) greater in among
  ## the "treated" (however defined).
  ##
  ## It takes five arguments:
  ##   'totalN': total number of observations
  ##   'treatedN': number of observations that received treatment
  ##   'totalSuccesses': total number of "successes"--i.e., 1's
  ##   'treatedSuccesses': number of successes among the treated
  ##   'Gammas': a vector Gammas (bounds on the differential odds of 
  ##  treatment) at which to test the significance of the results.
  ##
  ## It returns a matrix of Gammas and upper and lower bounds on the exact 
  ## p-value for Fisher's test.
  
  n <- totalN
  m <- treatedN
  c_plus <- totalSuccesses
  a <- treatedSuccesses
  
  Upsilon <- function(n, m, c_plus, a, Gamma) {
    ## Prob A >= a 
    numer <- 0
    for(k in max(a, (m + c_plus - n)):min(m, c_plus)) {
      numer <- numer +
        choose(c_plus, k) * choose((n-c_plus), m - k) * (Gamma ^ k)
    }
    denom <- 0
    for(k in max(0, (m + c_plus - n)):min(m, c_plus)) {
      denom <- denom +
        choose(c_plus, k) * choose((n - c_plus), (m - k)) * Gamma ^ k
    }
    return(numer / denom)
  }
  p_plus <- rep(NA, length(Gammas))
  p_minus <- rep(NA, length(Gammas))
  
  for(g in 1:length(Gammas)) {
    p_plus[g] <- Upsilon(n = n, m = m, c_plus = c_plus,
                         a = a, Gamma = Gammas[g])
    p_minus[g] <- Upsilon(n = n, m = m, c_plus = c_plus,
                          a = a, Gamma = (1 / Gammas[g]))
  }
  
  output <- cbind(Gammas, p_minus, p_plus)
  dimnames(output)[[2]] <- c("Gamma", "P-Value LB", "P-Value UB")
  return(output)
}
