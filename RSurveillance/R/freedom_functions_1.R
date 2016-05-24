##############################################################################################
# Freedom functions
##############################################################################################

# Population sensitivity calculations
###########################################################################
# sep.binom
# sep.hypergeo
# sep.exact
# spp
# sep
# sep.var.se
# sep.sys

##' Binomial Population sensitivity
##' @description Calculates population sensitivity for detecting disease,
##'   assuming imperfect test sensitivity and specificity and representative sampling,
##'   using binomial distribution (assumes large or unknown population size and that 
##'   cut-point number of reactors for a positive result = 1) 
##' @param n sample size = number of units tested (integer), scalar or vector  
##' @param pstar design prevalence as a proportion (scalar or vector of same length as n)
##' @param se unit sensitivity of test (proportion), default = 1 (scalar or vector of same length as n)
##' @param sp unit specificity of test (proportion), default = 1 (scalar or vector of same length as n)
##' @return vector of population-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.binom - checked
##' sep.binom(n=300, pstar = 0.02, se = 0.92)
##' tested<- seq(10,100, by=10)
##' prev<- 0.05
##' sens<- 0.9
##' sep.binom(tested, prev, sens)
sep.binom<- function(n, pstar, se = 1, sp = 1) {
    sep<- 1 - (1 - (se*pstar+(1-sp)*(1-pstar)))^n
    return(sep)
} # end of sep.binom

##' Population sensitivity for census (all units tested)
##' @description Calculates population sensitivity for detecting disease
##'   assuming imperfect test sensitivity, perfect test specificity
##'   and a census of all units in the population
##' @param se unit sensitivity of test (proportion), scalar or vector
##' @param d expected number of infected units in population (=design prevalence*N
##'   rounded to next integer), scalar or vector of same length as se
##' @return vector of population-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.exact - checked
##' sep.exact(d=1, se = 0.92)
##' inf<- 1:5
##' sens<- 0.8
##' sep.exact(d=inf, se=sens)
##' sep.exact(se=0.8, d = ceiling(0.01*c(10, 50, 100, 250, 500)))
sep.exact<- function(d=1, se = 1) {
  sep<- 1 - (1 - se)^d
    return(sep)
} # end of sep.exact


##' Hypergeometric Population sensitivity
##' @description Calculates population sensitivity for detecting disease,
##'   assuming imperfect test sensitivity, perfect test specificity 
##'   and representative sampling,
##'   using hypergeometric approximation (assumes known population size) 
##' @param N population size, scalar or vector of same length as n
##' @param n sample size (number tested), scalar or vector
##' @param d expected number of infected units in population (=design prevalence*N 
##'   rounded to next integer)
##' @param se unit sensitivity of test (proportion), scalar or vector of same length as n
##' @return a vector of population-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.hypergeo - checked
##' sep.hypergeo(N=100, n=50, d=1, se = 0.92)
##' inf<- 1:5
##' sens<- 0.8
##' sep.hypergeo(N=100, n=50, d=inf, se=sens)
##' N<- c(10, 50, 100, 250, 500)
##' sep.hypergeo(se=0.8, N=N, n=c(5, 25, 50, 125, 250), d = ceiling(0.01*N))
sep.hypergeo<- function(N, n, d, se = 1) {
  sep<- 1 - (1 - se*n/N)^d
    return(sep)
} # end of sep.hypergeo


##' Population specificity
##' @description Calculates population specificity assuming representative sampling
##' @param n sample size (number tested), integer, scalar or vector  
##' @param sp unit specificity of test (proportion), scalar or vector of same length as n
##' @return a vector of population-level specificities
##' @keywords methods
##' @export
##' @examples 
##' # examples for spp - checked
##' spp(10, 0.9)
##' spp(c(10, 20, 50, 100), 0.99)
##' spp(100, c(0.999, 0.99, 0.98, 0.95, 0.9))
spp<- function(n, sp) {
  sph<- sp^n
  return(sph)
} # end of spp function


##' Population sensitivity 
##' @description Calculates population sensitivity using appropriate method,
##'   depending on whether or not N provided (hypergeometric if N provided, 
##'   binomial otherwise), assuming perfect 
##'   test specificity and representative sampling
##' @param N population size, NA or vector of same length as n
##' @param n sample size (number tested), scalar or vector
##' @param pstar design prevalence as a proportion or integer, scalar
##'   or vector of same length as n
##' @param se unit sensitivity, scalar or vector of same length as n
##' @return a vector of population-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep - checked
##' sep(n=300, pstar=0.01, se=1)
##' sep(NA, 300, 0.01, 1)
##' sep(10000, 150, 0.02, 1)
##' sep(n=1:100, pstar = 0.05, se=0.95)
##' N<- seq(30, 100, by = 5)
##' se<- 0.95
##' pstar<- 0.1
##' n<- rep(30, length(N))
##' sep(N, n, pstar, se = se)
##' sep(rep(100, 10), seq(10, 100, by = 10), pstar = 1, se=0.99)
##' N<- c(55, 134, NA, 44, 256)
##' n<- c(15, 30, 28, 15, 33)
##' sep(N, n, 0.1, 0.95)
sep<- function(N = NA, n, pstar, se=1) {
  # check for errors in inputs
  # pstar.int = flag to indicate proportion (F) or integer (T) design prevalence
  pstar.int<- !(pstar < 1 & pstar > 0)
  if (sum(is.na(N)) > 0 & pstar.int) {
    err.msg<- "Population size (N) must be provided if design prevalence is an integer."
    return(err.msg)
  } else if (sum(n[!is.na(N)] > N[!is.na(N)]) > 0) {
    err.msg<- "Sample size (n) cannot be greater than population size (N)."
    return(err.msg)
  } else if (pstar.int & (pstar < 1 | pstar != round(pstar, 0))) {
    err.msg<- "Design prevalence must be an integer >= 1 if it is specified as an integer."
    return(err.msg)
  } else if (!pstar.int & (pstar >= 1 | pstar <= 0)) {
    err.msg<- "Design prevalence must be >0 and <1 if it is specified as a proportion."
    return(err.msg)
  }
  # sep calculations
  sep<- numeric(length(n))
  if (length(N) == 1) N<- rep(N, length(n))
  if (length(se) == 1) se<- rep(se, length(n))
  d<- pstar
  if (length(d) == 1) d<- rep(d, length(n))
  if (length(sep[is.na(N)]) >0) sep[is.na(N)]<- sep.binom(n=n[is.na(N)], pstar=pstar, se=se[is.na(N)])
  if (sum(!is.na(N)) != 0) {
    if (!pstar.int) {
      d[!is.na(N)]<- ceiling(N[!is.na(N)] * pstar)
    }
    sep[N == n & !is.na(N)]<- sep.exact(d=d[N == n & !is.na(N)], se=se[N == n & !is.na(N)])
    sep[N != n & !is.na(N)]<- sep.hypergeo(N=N[N != n & !is.na(N)], n=n[N != n & !is.na(N)], d=d[N != n & !is.na(N)], se=se[N != n & !is.na(N)])

  }
  return(sep)
}


##' Population sensitivity for varying unit sensitivity
##' @description Calculates population-level sensitivity where unit sensitivity 
##'   varies and using the appropriate method, depending on whether or not N provided 
##'   (hypergeometric if N provided, binomial otherwise), assuming perfect 
##'   test specificity and representative sampling
##' @param N population size (number of units or clusters), N must be >= length(se)) 
##'   or NA if unknown
##' @param se vector of unit sensitivity values (proportion) for each unit sampled
##' @param pstar specified design prevalence (scalar)
##' @return a scalar of population-level sensitivity
##' @keywords methods
##' @export
##' @examples 
##' # examples of sep.var.se - checked
##' sens<- c(rep(0.9, 50), rep(0.95, 100))
##' sep.var.se(NA, sens, 0.01)
##' sep.var.se(se=sens, pstar=0.01)
##' sep.var.se(N=500, sens, 0.01)
##' sep.var.se(NA, runif(150, 0.95, 0.99), 0.02)
##' sep.var.se(500, runif(150, 0.95, 0.99), 0.02)
sep.var.se<- function(N=NA, se, pstar) {
  if (!(is.na(N)) & N < length(se)) return("Error: N cannot be less than the number of sensitivity values")
  if (is.na(N)) {
    sep<- 1-prod(1-se*pstar)
  } else {
    if (pstar < 1 & pstar > 0) pstar<- ceiling(N*pstar)
    sep<- 1-(1-mean(se)*length(se)/N)^pstar
  }
  return(sep)
}


##' 2-stage population sensitivity 
##' @description Calculates population-level (system) sensitivity for representative
##'   2-stage sampling (sampling of clusters and units within clusters), 
##'   assuming imperfect test sensitivity and perfect test specificity 
##' @param H population size = number of clusters in the population, default = NA
##' @param N population size within clusters,
##'   scalar or a vector of same length as n, default = NA
##' @param n sample size (vector of number tested per cluster)
##' @param se unit sensitivity of test (proportion), scalar, default = 1
##' @param pstar.c cluster (herd) level design prevalence, scalar,
##'   either proportion or integer
##' @param pstar.u unit (animal) level design prevalence, scalar,
##'   either proportion or integer
##' @return list of 6 elements, 1) population level sensitivity, 2) vector of 
##'   cluster-level sensitivities, 3) N, 4) n, 5) vector of design prevalences 
##'   and 6) unit sensitivity 
##' @note if pstar.c is not a proportion N must be provided 
##'   (and N>=n)
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.sys - checked
##' H<- 500
##' N<- rep(1000, 150)
##' N[5]<- NA
##' n<- rep(30, 150)
##' pstar.u<- 0.1
##' pstar.c<- 0.01
##' se<- 0.98
##' sep.sys(H, N, n, pstar.c, pstar.u, se)
##' sep.sys(NA, N, n, 0.02, 0.05, 0.95)
##' N<- round(runif(105)*900+100)
##' n<- round(runif(105)*30+10)
##' sse<- sep.sys(1000, N, n, 0.02, 0.05, 0.9)
##' data.frame(N, n, sse[[2]])
sep.sys<- function(H=NA, N=NA, n, pstar.c, pstar.u, se=1) {
  # calculate cluster level sensitivities
  sep.cluster<- sep(N, n, pstar.u, se)
  # calculate overall system sensitivity
  sep<- sep.var.se(H, sep.cluster, pstar.c)
  return(list("Population_sensitivity"=sep,
              "Cluster_sensitivities"=sep.cluster,
              "N" = N,
              "n"= n,
              "Design prevalence"= c("cluster-level"=pstar.c, "unit-level"=pstar.u),
              "Unit sensitivity"= se))
}


###########################################################################
# Sample size calculations
###########################################################################
# calculate sample size assuming sampling with (binomial) or 
# without (hypergeometric) replacement
# n.binom
# n.hypergeo
# n.freedom
# n.2stage

##' Binomial sample size
##'@description Calculates sample size for demonstrating freedom or 
##'   detecting disease using binomial approach and assuming 
##'   imperfect test sensitivity, perfect test specificity and
##'   representative sampling
##' @param sep desired population sensitivity (scalar or vector)
##' @param pstar specified design prevalence (scalar or vector of same length as sep)
##' @param se unit sensitivity, default = 1 (scalar or vector of same length as sep)
##' @return vector of sample sizes
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.binom - checked
##' n.binom(sep=0.95, pstar=c(0.01, 0.02, 0.05, 0.1, 0.2))
##' n.binom(c(0.5, 0.8, 0.9, 0.95), 0.01)
n.binom<- function(sep, pstar, se = 1) {
  n<- log(1 - sep)/log(1 - pstar * se)
  return(ceiling(n))
}


##' Hypergeometric sample size
##'@description Calculates sample size for demonstrating freedom or 
##'   detecting disease using hypergeometric approximation and assuming 
##'   imperfect test sensitivity, perfect test specificity and
##'   representative sampling
##' @param sep desired population sensitivity (scalar or vector)
##' @param N population size (scalar or vector of same length as sep)
##' @param d expected number of infected units in population, = design prevalence*N
##'   rounded to next integer (scalar or vector of same length as sep)
##' @param se unit sensitivity, default = 1 (scalar or vector of same length as sep)
##' @return vector of sample sizes, NA if n>N
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.hypergeo - checked
##' n.hypergeo(0.95, N=100, d=1, se = 0.95)
##' n.hypergeo(sep=0.95, N=c(100, 200, 500, 1000, 10000), d=ceiling(0.01*c(100, 200, 500, 1000, 10000)))
##' n.hypergeo(c(0.5, 0.8, 0.9, 0.95), N=100, d=5)
##' n.hypergeo(0.95, N=80, d=c(1, 2, 5, 10))
##' n.hypergeo(0.95, N=80, d=c(1, 2, 5, 10), se = 0.8)
n.hypergeo<- function(sep, N, d, se = 1) {
  n<- (N/se)*(1 - (1 - sep)^(1/d))
  n[n > N]<- NA
  return(ceiling(n))
}


##' Freedom sample size 
##' @description Calculates sample size for demonstrating freedom or 
##'   detecting disease using the appropriate method, depending on 
##'   whether or not N provided (hypergeometric if N provided, binomial otherwise), 
##'   assuming imperfect test sensitivity, perfect test specificity 
##'   and representative sampling
##' @param N population size, default = NA (unknown) (scalar or vector of same length as sep)
##' @param sep desired population sensitivity (scalar or vector)
##' @param pstar specified design prevalence as proportion or integer 
##'   (scalar or vector of same length as sep)
##' @param se unit sensitivity (scalar or vector of same length as sep)
##' @return vector of sample sizes, NA if N is specified and n>N
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.freedom - checked
##' n.freedom(NA, sep=0.95, pstar=0.01, se=1)
##' n.freedom(500, sep=0.95, pstar=0.01, se=1)
##' n.freedom(N=c(100, 500, 1000, 5000, 10000, 100000, NA), sep=0.95, pstar=0.01, se=1)
##' n.freedom(500, sep=0.95, pstar=0.01, se=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1))
n.freedom<- function(N=NA, sep=0.95, pstar,se=1) {
  if (length(N) == 1) {
    if (is.na(N)) {
      n<- n.binom(sep, pstar, se)
    } else {
      d<- pstar
      if (pstar < 1 & pstar > 0) {
        d<- ceiling(N * pstar)
      }
      n<- n.hypergeo(sep, N, d, se)
    }
  } else {
    n<- numeric(length(N))
    n[is.na(N)]<- n.binom(sep=sep, pstar=pstar, se=se)
    pstar.int<- !(pstar < 1 & pstar > 0)
    d<- pstar
    if (length(d) == 1) d<- rep(d, length(N))
    if (pstar < 1 & pstar > 0) {
      d[!is.na(N)]<- ceiling(N[!is.na(N)] * pstar)
    }
    n[!is.na(N)]<- n.hypergeo(sep, N[!is.na(N)], d[!is.na(N)], se)
  }
  return(n)
}


##' 2-stage freedom sample size
##' @description Calculates sample sizes for a 2-stage representative survey 
##'   (sampling of clusters and units within clusters) for disease freedom or detection,
##'   assuming imperfect test sensitivity, perfect test specificity and representative sampling
##' @param H population size = number of clusters or NA if not known, default = NA
##' @param N population sizes for clusters, default = NA, scalar or 
##'   vector of population sizes for clusters
##' @param sep.sys desired population sensitivity (scalar)
##' @param sep.c desired cluster-level sensitivity (scalar)
##' @param pstar.c specified cluster-level design prevalence as 
##'   proportion or integer (scalar) 
##' @param pstar.u specified population-level design prevalence as 
##'   proportion or integer (scalar)
##' @param se unit sensitivity (scalar)
##' @return a list of 2 elements, the number of clusters to sample and a vector of 
##'   sample sizes per cluster
##' @keywords methods
##' @export
##' @examples 
##' # examples of n.2stage - checked
##' n.2stage(NA, NA, 0.95, 0.5, 0.01, 0.1, 0.95)
##' n.2stage(500, NA, 0.95, 0.5, 10, 0.1, 0.95)
##' n.2stage(1000, c(50, 100, 200, 500, 1000, 5000, NA), 0.95, 0.5, 0.01, 0.05, 0.8)
##' n.2stage(1000, c(50, 100, 200, 500, 1000, 5000, NA), 0.95, 0.5, 0.01, 1, 0.8)
##' n.2stage(1000, c(50, 100, 200, 500, 1000, 5000, NA), 0.9, 0.95, 1, 0.1, 0.8)
n.2stage<- function(H=NA, N=NA, sep.sys=0.95, sep.c, pstar.c, pstar.u, se=1) {
  n.clusters<- n.freedom(H, sep.sys, pstar.c, sep.c)
  n.units<- n.freedom(N, sep.c, pstar.u, se)
  names(n.units)<- N
  names(n.units)[length(N)]<- "Inf"
  return(list("Clusters" = n.clusters, "Units" = n.units))
}


###########################################################################
# design prevalence for given sample size
###########################################################################
# pstar.calc

##' Design prevalence back-calculation
##' @description Calculates design prevalence required for given sample size and
##'   desired population-level sensitivity, assuming 
##'   imperfect test sensitivity, perfect test specificity and
##'   representative sampling
##' @param N populaton size if known (scalar or vector of same length as n)
##' @param n sample size (scalar or vector)
##' @param sep desired population sensitivity (scalar or vector of same length as n)
##' @param se unit sensitivity (scalar or vector of same length as n)
##' @return vector of design prevalence values
##' @keywords methods
##' @export
##' @examples 
##' # examples of pstar.calc- checked
##' pstar.calc(NA, 280, 0.95, 0.98)
##' pstar.calc(500, 250, sep=0.95, se=1)
##' pstar.calc(N=c(100, 500, 1000, 5000, 10000, 100000, NA), n=30, sep=0.95, se=1)
##' pstar.calc(500, n=30, sep=0.95, se=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1))
pstar.calc<- function(N=NA, n, sep, se) {
  if (length(N) == 1) {
    if (is.na(N)) {
      pstar<- (1 - exp((log(1 - sep))/n))/se
    } else {
      pstar<- log(1 - sep)/log(1 - se * n/N)/N
    }
  } else {
    if (length(n) == 1) n<- rep(n, length(N))
    pstar<- numeric(length(N))
    pstar[is.na(N)]<- (1 - exp((log(1 - sep))/n[is.na(N)]))/se
    pstar[!is.na(N)]<- log(1 - sep)/log(1 - se * n[!is.na(N)]/N[!is.na(N)])/N[!is.na(N)]
  }
  return(pstar)
}


###########################################################################

# Confidence of freedom
###########################################################################
# pfree.equ
# disc.prior
# pfree.1
# pfree.calc
# sep.pfree
# sep.prior
# n.pfree

##' Equilibrium probability of freedom
##' @description Calculates equilibrium probability of disease freedom and 
##'   equilibrium prior probability of freedom, after discounting for
##'   probability of introduction
##' @param sep population sensitivity for time period (scalar or
##'   vector)
##' @param p.intro probability of introduction for time period (scalar
##'   or vector of same length as sep)
##' @return a list of 2 vectors, equilibrium posterior probability of freedom 
##'   and equilibrium prior (discounted) probability of freedom
##' @keywords methods
##' @export
##' @examples 
##' # examples of pfree.equ
##' pfree.equ(runif(10, 0.4, 0.6), 0.01)
##' pfree.equ(0.8, 0.05)
##' pfree.equ(rep(0.9, 6), c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05))
pfree.equ<- function(sep, p.intro) {
  pf.equ<- (1 - (p.intro/sep))/(1 - p.intro)
  prior.equ<- 1 - (p.intro/sep)
  return(list("Equ_PFree"= pf.equ, "Equ_prior"= prior.equ))
}


##' Discounted prior probability of freedom
##' @description Calculates the discounted prior probability of disease freedom,
##'   after adjusting for the probability of disease exceeding the 
##'   design prevalence during the time period of the surveillance data being analysed
##' @param prior prior probability of freedom before surveillance
##' @param p.intro probability of introduction 
##'   (or of prevalence exceeding the design prevalence) during the time period 
##'   (scalar or vector equal length to prior)
##' @return vector of discounted prior probabilities of freedom
##' @keywords methods
##' @export
##' @examples 
##' # examples for disc.prior
##' disc.prior(0.5, 0.01)
##' disc.prior(0.95, c(0.001, 0.005, 0.01, 0.02, 0.05))
##' disc.prior(c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95), 0.01)
disc.prior<- function(prior, p.intro) {
  prior.disc<- 1 - (1 - prior + p.intro - ((1 - prior) * p.intro))
  return(prior.disc)
}


##' Probability of freedom for single time period
##' @description Calculates the posterior probability (confidence) of disease 
##'   freedom (negative predictive value) for a single time period
##' @param sep population sensitivity for time period (scalar or
##'   vector)
##' @param p.intro probability of introduction for time period (scalar
##'   or vector of same length as sep)
##' @param prior prior probability of freedom before surveillance (scalar
##'   or vector of same length as sep)
##' @return \code{data.frame} with columns for sep, p.intro, discounted
##' prior, pfree, pfree.equ and prior.equ
##' @keywords methods
##' @export
##' @examples 
##' # examples for pfree.1
##' pfree.1(0.8, 0.01, 0.5)
##' pfree.1(0.6, c(0.001, 0.005, 0.01, 0.02, 0.05), 0.5)
##' pfree.1(runif(10, 0.4, 0.6), 0.01, 0.5)
##' pfree.1(runif(10, 0.4, 0.6), runif(10, 0.005, 0.015), 0.5)
pfree.1<- function(sep, p.intro, prior=0.5) {
  if (length(p.intro) < length(sep)) p.intro<- rep(p.intro, length(sep))
  prior.disc<- numeric(length(sep))
  pfree<- numeric(length(sep))
  prior.disc<- disc.prior(prior, p.intro)
  pfree<- prior.disc/(1 - sep * (1-prior.disc))
  tmp<- pfree.equ(sep, p.intro)
  pfree.eq<- tmp[[1]]
  prior.eq<- tmp[[2]]
  return(data.frame(SeP=sep,
                    PIntro=p.intro,
                    "Discounted prior"=prior.disc,
                    PFree=pfree,
                    "Equilibrium PFree"= pfree.eq,
                    "Equilibrium prior"=prior.eq))
}


##' Probability of freedom over time
##' @description Calculates the probability (confidence) of disease freedom for 
##' given prior, sep and p.intro over 1 or more time periods
##' @param sep population sensitivity for each time period (vector)
##' @param p.intro probability of introduction for each time period (scalar
##'   or vector of same length as sep)
##' @param prior prior probability of freedom before surveillance (scalar)
##' @return \code{data.frame} with columns for sep, p.intro, discounted
##'   prior, probability of freedom, equilibrium probability of freedom
##'   and equilibrium prior
##' @keywords methods
##' @export
##' @examples 
##' # examples for pfree.calc
##' pfree.calc(0.8, 0.01, 0.5)
##' pfree.calc(rep(0.6,24), 0.01, 0.5)
##' pfree.calc(runif(10, 0.4, 0.6), 0.01, 0.5)
##' pfree.calc(runif(10, 0.4, 0.6), runif(10, 0.005, 0.015), 0.5)
pfree.calc<- function(sep, p.intro, prior=0.5) {
  if (length(p.intro) < length(sep)) p.intro<- rep(p.intro, length(sep))
  prior.disc<- numeric(length(sep))
  pfree<- numeric(length(sep))
  pfree[1]<- pfree.1(sep[1], p.intro[1], prior[1])[,4]
  prior.disc[1]<- disc.prior(prior, p.intro[1])
  if (length(sep) > 1) {
    for (p in 2:length(sep)) {
      prior.disc[p]<- disc.prior(pfree[p-1], p.intro[p])
      pfree[p]<- pfree.1(sep[p], p.intro[p], pfree[p-1])[,4]
    }
  }
  tmp<- pfree.equ(sep, p.intro)
  pfree.eq<- tmp[[1]]
  prior.eq<- tmp[[2]]
  return(data.frame("Period"=1:length(sep),
                    SeP=sep,
                    PIntro=p.intro,
                    "Discounted prior"=prior.disc,
                    PFree=pfree,
                    "Equilibrium PFree"= pfree.eq,
                    "Equilibrium prior"=prior.eq))
}


##' Population sensitivity to achieve desired (posterior) probability of freedom
##' @description Calculates the population sensitivity required to achieve a 
##'   given value for probability of disease freedom
##' @param prior prior probability of freedom before surveillance (scalar or vector)
##' @param pfree desired probability of freedom (scalar or vector)
##' @return a vector of population-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples of sep.pfree
##' sep.pfree(0.5, 0.95)
##' sep.pfree(c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95), 0.99)
##' sep.pfree(0.5, c(0.8, 0.9, 0.95, 0.99))
sep.pfree<- function(prior, pfree) {
  sep<- (1-prior/pfree)/(1-prior)
  return (sep)
}


##' Population sensitivity to achieve desired prior probability of freedom
##' @description Calculates the population sensitivity required to achieve a 
##'   given value for the prior (discounted) probability of disease freedom
##' @param prior prior probability of freedom before surveillance (scalar or vector)
##' @param p.intro probability of introduction for time period (scalar or vector equal length to sep)
##' @return a vector of population-level sensitivities
##' @keywords methods
##' @export
##' @examples 
##' # examples of sep.prior
##' sep.prior(0.95, 0.01)
##' sep.prior(c(0.9, 0.95, 0.98, 0.99), 0.01)
##' sep.prior(0.95, c(0.001, 0.005, 0.01, 0.02, 0.05))
sep.prior<- function(prior, p.intro) {
  sep<- p.intro/(1-prior)
  return (sep)
}


##' Sample size to achieve desired (posterior) probability of freedom
##' @description Calculates the sample size required to achieve a 
##'   given value for probability of disease freedom
##' @param pfree desired probability of freedom (scalar or vector)
##' @param prior prior probability of freedom before surveillance (scalar
##'   or vector of same length as pfree)
##' @param p.intro probability of introduction for time period (scalar
##'   or vector of same length as pfree)
##' @param se unit sensitivity (scalar
##' or vector of same length as pfree)
##' @param pstar design prevalence (scalar
##' or vector of same length as pfree)
##' @param N population size (scalar
##' or vector of same length as pfree)
##' @return vector of sample sizes
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.pfree
##' n.pfree(0.95, 0.5, 0.01, 0.05, 0.9)
##' n.pfree(0.95, 0.5, 0.01, 0.05, 0.9, N=300)
##' n.pfree(pfree = c(0.9, 0.95, 0.98, 0.99), prior = 0.7, 0.01, 0.01, 0.8, 1000)
##' n.pfree(0.95, 0.7, 0.01, 0.1, 0.96)
n.pfree<- function(pfree, prior, p.intro, pstar, se, N = NA) {
  sep<- sep.pfree(prior, pfree)
  n<- n.freedom(N, sep, pstar, se)
  return(n)
}


###########################################################################
# combining tests
###########################################################################
# se.series
# se.parallel
# sp.series
# sp.parallel

##' Sensitivity of tests in series
##' @description Calculates the combined sensitivity for multiple tests 
##'   interpreted in series (assuming independence)
##' @param se vector of unit sensitivity values
##' @return scalar of combined sensitivity, assuming independence
##' @keywords methods
##' @export
##' @examples 
##' # examples for se.series
##' se.series(c(0.99, 0.95, 0.8))
se.series<- function(se) {
  se.comb<- prod(se)
  return(se.comb)
}


##' Sensitivity of tests in parallel
##' @description Calculates the combined sensitivity for multiple tests 
##'   interpreted in parallel (assuming independence)
##' @param se vector of unit sensitivity values
##' @return scalar of combined sensitivity, assuming independence
##' @keywords methods
##' @export
##' @examples 
##' # examples for se.parallel
##' se.parallel(c(0.99, 0.95, 0.8))
se.parallel<- function(se) {
  se.comb<- 1 - prod(1-se)
  return(se.comb)
}

##' Specficity of tests in series
##' @description Calculates the combined specificity for multiple tests 
##'   interpreted in series (assuming independence)
##' @param sp vector of unit specificity values
##' @return scalar of combined specificity, assuming independence
##' @keywords methods
##' @export
##' @examples 
##' # examples for sp.series
##' sp.series(c(0.99, 0.95, 0.8))
sp.series<- function(sp) {
  sp.comb<- 1 - prod(1 - sp)
  return(sp.comb)
}


##' Specificity of tests in parallel
##' @description Calculates the combined specificity for multiple tests 
##'   interpreted in parallel (assuming independence)
##' @param sp vector of unit specificity values
##' @return scalar of combined specificity, assuming independence
##' @keywords methods
##' @export
##' @examples 
##' # examples for sp.parallel
##' sp.parallel(c(0.99, 0.95, 0.8))
sp.parallel<- function(sp) {
  sp.comb<- prod(sp)
  return(sp.comb)
}

############################################################################
# pooled sampling (see Christensen & Gardner (2000). Herd-level interpretation of test results for
# epidemiologic studies of animal diseases. Prev Vet Med. <b>45</b>:83-106
############################################################################
# sep.pooled
# n.pooled

##' Pooled population sensitivity
##' @description Calculates population sensitivity (sep) and population specificity (spp)
##'   assuming pooled sampling
##'   and allowing for imperfect sensitivity and specificity of the pooled test
##' @param r number of pools sampled (scalar or vector)
##' @param k pool size (scalar or vector of same length as r)
##' @param pstar design prevalence (scalar or vector of same length as r)
##' @param pse pool-level sensitivity (scalar or vector of same length as r)
##' @param psp pool-level specificity (scalar or vector of same length as r)
##' @return list of 2 elements, vector of sep values and vector of spp
##'   values
##' @keywords methods
##' @export
##' @examples 
##' # examples for sep.pooled
##' sep.pooled(60, 5, 0.01, 1, 1)
##' sep.pooled(4, 10, 0.1, 0.9, 1)
##' sep.pooled(1:10*5, 5, 0.02, 0.9, 0.99)
##' sep.pooled(10, 5, 0.05, c(0.8, 0.9, 0.95, 0.99), 1)
sep.pooled<- function(r, k, pstar, pse, psp=1) {
  sep<- 1 - ((1 - (1 - pstar)^k)*(1 - pse) + (1 - pstar)^k * psp)^r
  spp<- psp^r
  return(list(sep=sep, spp=spp))
}


##' Sample size for pooled testing for freedom
##' @description Calculates sample size to achieve desired 
##'   population-level sensitivity, assuming pooled sampling
##'   and allowing for imperfect sensitivity and specificity of the pooled test
##' @param sep desired population sensitivity (scalar or vector)
##' @param k pool size (constant across pools) (scalar or vector of same length as sep)
##' @param pstar design prevalence (scalar or vector of same length as sep)
##' @param pse pool-level sensitivity (scalar or vector of same length as sep)
##' @param psp pool-level specificity (scalar or vector of same length as sep)
##' @return vector of sample sizes
##' @keywords methods
##' @export
##' @examples 
##' # examples for n.pooled
##' n.pooled(0.95, 5, 0.01, 1, 1)
##' n.pooled(0.95, 10, 0.1, 0.9, 1)
##' n.pooled(0.95, c(2, 5, 10, 20), 0.1, c(0.99, 0.98, 0.97, 0.95), 1)
n.pooled<- function(sep, k, pstar, pse, psp=1) {
  n<- log(1-sep)/log(((1 - (1 - pstar)^k)*(1 - pse) + (1 - pstar)^k * psp))
  return(ceiling(n))
}
