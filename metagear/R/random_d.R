#' Random generation of Hedges' d effect sizes.  
#'
#' Generates random Hedges' d (1981, 1982) effect sizes and their variances.       
#'
#' @param K Number of effect sizes to generate.
#' @param X_t The population mean (mu) of the (t)reatment group.
#' @param var_t The population variance of the treatment group mean.
#' @param N_t The number of samples of the treatment mean.  When a non-negative
#'    integer, all treatment means will be estimated using the same N.  A
#'    vector of unequal N's can also be taken; if so, K will be ignored and the 
#'    number of randomly generated means will equal the length of that vector,
#'    and each mean will be based on each N within the vector.  
#' @param X_c The population mean (mu) of the (c)ontrol group.
#' @param var_c The population variance of the control group mean.
#' @param N_c The number of samples of the control mean.  When a non-negative
#'    integer, all control means will be estimated using the same N.  A
#'    vector of unequal N's can also be taken; if so, K will be ignored and the 
#'    number of randomly generated means will equal the length of that vector,
#'    and each mean will be based on each N within the vector.  
#' @param bias_correction When \code{"FALSE"}, returns Cohen's g effect sizes 
#'    that are not adjusted using a small-sample correction (J).
#'
#' @return A data table with columns of random effect sizes (d) and their 
#'    variances (var_d).
#'
#' @examples
#'    random_d(K = 5, X_t = 25, var_t = 1, N_t = 15, X_c = 10, var_c = 1, N_c = 15)
#'
#' @references Hedges, L.V. 1981. Distribution theory for Glass's estimator of 
#'    effect size and related estimators. Journal of Educational Statistics 
#'    6: 107-128.
#' @references Hedges, L.V. 1982. Estimation of effect size from a series of 
#'    independent experiments. Psychological Bulletin 92: 490-499.
#'
#' @importFrom stats rnorm sd
#' @export random_d

random_d <- function(K,
                     X_t,
                     var_t,
                     N_t,
                     X_c,
                     var_c,
                     N_c,
                     bias_correction = TRUE) {
  
  # add error if N_t & N_c are not same length
  
  if(length(N_t) != 1) {
    samples_t <- sum(N_t)
    theGroups_t <- rep(1:length(N_t), N_t)
    #add warning message that K is ignored
  }
  else {
    samples_t <- K * N_t
    theGroups_t <- rep(1:K, each = N_t)
  }
  
  if(length(N_c) != 1) {
    samples_c <- sum(N_c)
    theGroups_c <- rep(1:length(N_c), N_c)
    #add warning message that K is ignored
  }
  else {
    samples_c <- K * N_c
    theGroups_c <- rep(1:K, each = N_c)
  }
  
  tList <- data.frame(X = rnorm(samples_t, 
                                mean = X_t, 
                                sd = sqrt(var_t)),
                      group = theGroups_t)  
  cList <- data.frame(Y = rnorm(samples_c, 
                                mean = X_c, 
                                sd = sqrt(var_c)), 
                      group = theGroups_c)  
  
  means_t <- tapply(tList$X, list(tList$group), mean)
  sds_t <- tapply(tList$X, list(tList$group), sd)
  means_c <- tapply(cList$Y, list(cList$group), mean)
  sds_c <- tapply(cList$Y, list(cList$group), sd);
  
  # add non approximation form
  if(bias_correction == TRUE) {
    J <- 1 - 3 / (4 * (N_c + N_t - 2) - 1)
  } else {
    J <- 1.0
  }
  
  pooled_SD <- sqrt(((N_c - 1)* sds_c ^ 2 + (N_t - 1)* sds_t ^ 2) / (N_c + N_t - 2))
  d <- ((means_t - means_c) / pooled_SD) * J
  
  thed <- data.frame(d = d, 
                     var_d = (N_c + N_t) / (N_c * N_t) + (d ^ 2) / (2 * (N_c + N_t)))
  
  return(thed)  
}
