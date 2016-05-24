#' Random generation of log response ratio (RR) effect sizes.  
#'
#' Generates random log response ratios and their variances (Hedges et al. 1999).
#' NOTE: samples from a log-normal distribution to generate non-negative control
#' and treatment means (following Lajeunesse 2015).       
#'
#' @param K Number of effect sizes to generate.
#' @param X_t The population mean (mu) of the (t)reatment group (numerator of 
#'    ratio).
#' @param var_t The population variance of the treatment group mean.
#' @param N_t The number of samples of the treatment mean.  When a non-negative
#'    integer, all treatment means will be estimated using the same N.  A
#'    vector of unequal N's can also be taken; if so, K will be ignored and the 
#'    number of randomly generated means will equal the length of that vector,
#'    and each mean will be based on each N within the vector.  
#' @param X_c The population mean (mu) of the (c)ontrol group (denominator of 
#'    ratio).
#' @param var_c The population variance of the control group mean.
#' @param N_c The number of samples of the control mean.  When a non-negative
#'    integer, all control means will be estimated using the same N.  A
#'    vector of unequal N's can also be taken; if so, K will be ignored and the 
#'    number of randomly generated means will equal the length of that vector,
#'    and each mean will be based on each N within the vector.  
#'
#' @return A data table with columns of random effect sizes (RR) and their 
#'    variances.
#'
#' @examples
#'    random_RR(K = 5, X_t = 25, var_t = 1, N_t = 15, X_c = 10, var_c = 1, N_c = 15)
#'
#' @references Hedges, L.V., J. Gurevitch, and P.S. Curtis. 1999. The 
#'    meta-analysis of response ratios in experimental ecology. Ecology 80: 
#'    1150-1156.
#' @references Lajeunesse, M.J. 2015. Bias and correction for the log response
#'    ratio used in ecological meta-analysis. Ecology.
#'
#' @importFrom stats rlnorm sd
#' @export random_RR

random_RR <- function(K,
                      X_t,
                      var_t,
                      N_t,
                      X_c,
                      var_c,
                      N_c) {
  
  # add error if N_t & N_c are not same length
  
  if(length(N_t) != 1) {
    samples_t <- sum(N_t)
    theGroups_t <- rep(1:length(N_t), N_t)
    finalN_t <- N_t
    #add warning message that K is ignored
  }
  else {
    samples_t <- K * N_t
    theGroups_t <- rep(1:K, each = N_t)
    finalN_t <- N_t
  }
  
  if(length(N_c) != 1) {
    samples_c <- sum(N_c)
    theGroups_c <- rep(1:length(N_c), N_c)
    finalN_c <- N_c
    #add warning message that K is ignored
  }
  else {
    samples_c <- K * N_c
    theGroups_c <- rep(1:K, each = N_c)
    finalN_c <- N_c
  }
  
  tList <- data.frame(X = rlnorm(samples_t, 
                                 meanlog = log((X_t^2)/(sqrt(var_t + X_t^2))), 
                                 sdlog = sqrt(log(var_t/(X_t^2) + 1.0))),
                      group = theGroups_t)  
  cList <- data.frame(Y = rlnorm(samples_c, 
                                 meanlog = log((X_c^2)/(sqrt(var_c + X_c^2))), 
                                 sdlog = sqrt(log(var_c/(X_c^2) + 1.0))), 
                      group = theGroups_c)  
  
  means_t <- tapply(tList$X, list(tList$group), mean)
  sds_t <- tapply(tList$X, list(tList$group), sd)
  means_c <- tapply(cList$Y, list(cList$group), mean)
  sds_c <- tapply(cList$Y, list(cList$group), sd);
    
  theRR <- data.frame(RR = log(means_t/means_c), 
                      var_RR = (sds_t^2)/(N_t * (means_t^2)) + (sds_c^2)/(N_c * (means_c^2)))  
  
  return(theRR)  
}
