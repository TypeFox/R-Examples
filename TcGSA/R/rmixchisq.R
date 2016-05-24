#'Random Generation of Chi-square Mixtures
#'
#'\code{rmixchisq} is used to simulate a mixture of chi-square distributions
#'that corresponds to the null distribution of the Likelihood Ratio between 2
#'nested mixed models.
#'
#'The approximate null distribution of a likelihood ratio for 2 nested mixed
#'models, where both fixed and random effects are tested simultaneously, is a
#'very specific mixture of \eqn{\chi^2}{\chi^2} distributions [\cite{Self & Liang
#'(1987), Stram & Lee (1994) and Stram & Lee (1995)}].  It depends on both the
#'number of random effects and the number of fixed effects to be tested
#'simultaneously: 
#'\deqn{LRT_{H_0}\sim\sum_{k=q}^{q+r}{{r}\choose{k-q}}2^{-r}\chi^2_{(k)}}{LRT_H0~\sum k=q..q+r combination(r,k-q) 2^(-r) \chi^2 (k)}
#'
#'@param n number of observations.
#'@param s number of fixed effects to be tested.
#'@param q number of random effects to be tested.
#'
#'@return A vector of random independent observations of the chisquare mixture
#'identified by the values of \code{s} and \code{q}.
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{pval_simu}}
#'
#'@references Self, S. G. and Liang, K., 1987, Asymptotic properties of maximum
#'likelihood estimators and likelihood ratio tests under nonstandard
#'conditions, \emph{Journal of the American Statistical Association} 82:
#'605--610.
#'
#'Stram, D. O. and Lee, J. W., 1994, Variance components testing in the
#'longitudinal mixed effects model, \emph{Biometrics} 50: 1171--1177.
#'
#'Stram, D. O. and Lee, J. W., 1995, Corrections to "Variance components
#'testing in the longitudinal mixed effects model" by Stram, D. O. and Lee, J.
#'W.; 50: 1171--1177 (1994), \emph{Biometrics} 51: 1196.
#'
#'
#'@importFrom stats rchisq
#'
#'@export
#'
#'@examples
#'library(graphics)
#'library(stats)
#'
#'sample_mixt <- rmixchisq(n=1000, s=3, q=3)
#'plot(density(sample_mixt))
#'
#'
rmixchisq <-
function(n,s,q){
  if(q>0){
    mixprobs <- numeric(q+1)
    for(k in (s:(q+s))){
      mixprobs[k-s+1] <- choose(q,k-s)*2^(-q)
    }
    mix <- n*mixprobs
    #s <- s + q*(q-1)/2 #conservative corrections for the covariances terms
  }else{
    mix=n
  }
  
  sample <- NULL
  for(k in (s:(q+s))){
    sample <- c(sample, stats::rchisq(mix[k-s+1],df=k))
  }
  
  return(sample) 
}
