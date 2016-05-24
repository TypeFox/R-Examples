#' Posterior Summary of Mortality
#' 
#' Calculates and summarizes the posterior distribution of mortality count.
#' 
#' Assuming a Gamma(xi, lam) on the average daily mortality rate m, this model
#' treats the mortality M for the current period as Poisson-distributed with
#' mean m*I. The carcass count C will include "new" carcasses with a Bi(M,T)
#' distribution as well as "old" carcasses (if bt > 0). For derivation of 
#' resulting conditional pdf see Wolpert (2015).
#' 
#' This function calls \code{acme.post} but suppresses plotting.
#' 
#'
#'@param C Observed mortality count. Non-negative integer or vector.
#'@param Rstar ACME inverse-inflation factor R*, reported by acme.summary() 
#'  as "Rstar."
#'@param T The first term in recursive calculation of Rstar, from acme.summary.
#'@param gam Values for highest posterior density credible interval.
#'@param I Interval length, days.
#'@param Mmax Maximimum value for which posterior probability is calculated.
#'@param xi First parameter of gamma prior. Default is 1/2 for Objective prior.
#'@param lam Second parameter of gamma prior. Default is 0 for Objective prior.
#'  
#'@export 
#'@return \code{acme.table} returns a table which includes ACME
#' estimate (M_hat), posterior mean, and highest posterior credible intervals for probabilities 
#' as specified by the parameter gam.
#' 
#' @examples
#' acme.table(C=0:5,Rstar = 0.2496, T = 0.174)
#' 

acme.table <- function(C=0, Rstar=0.2496, T=0.1740, gam=c(0.5, 0.9), I=7,
                       Mmax = 200, xi=1/2, lam=0){
  for(i in 1:length(C)){
    C_tab_i <- acme.post(C=C[i], Rstar=Rstar, T=T, gam=gam, I=I,
                         Mmax = Mmax, xi=xi, lam=lam, plotit=FALSE);
    if(i==1){C_tab <- C_tab_i}
    else{C_tab <- rbind(C_tab,C_tab_i)}
  }
  return(C_tab)
}