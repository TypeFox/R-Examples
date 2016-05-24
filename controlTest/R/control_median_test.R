#' Median Comparison for Two-Sample Right-Censored Survival Data
#' 
#' @param t1 A vector of observed right-censored survival times for group 1 (Control).
#' @param c1 A vector of censoring indicators for group 1 (0 = alive, 1 = dead).
#' @param t2 A vector of observed right-censored survival times for group 2 (Treatment).
#' @param c2 A vector of censoring indicators for group 2 (0 = alive, 1 = dead).
#' @param R Number of replications for bootstrapping (Default = 1000).
#' @param seed Seed number for bootstrapping (Default =  1234).
#' @return A list containing the median survival times for both groups, Z-score, and two-sided p-value.
#' A plot of the two survival curves is also outputted. 
#' @details This function only compares the median between two samples. A general quantile version has not been developed yet. It is important to note the 
#' possiblilty that the median survival time may not be estimable in our bootstrap samples. In such cases 
#' the largest observed survival time will be considered as an estimate for the median survival time. Also,
#' if the median survival time for the control group is larger than the longest survival time for the
#' treatment group, then Q will be evaluated using the last observed survival time for the treatment group.
#' @examples
#' #Reference: Klein and Moeschberger (1997) Survival Analysis Techniques for Censored and truncated data, Springer.
#' #Data: Chapter 7.6 Example 7.9 (p. 211)
#' 
#'library(controlTest)
#'t1 <- c(1, 63, 105, 129, 182, 216, 250, 262, 301, 301, 
#'        342, 354, 356, 358, 380, 383, 383, 338, 394, 408, 460, 489, 
#'        499, 523, 524, 535, 562, 569, 675, 676, 748, 778, 786, 797, 955, 968, 1000,
#'        1245, 1271, 1420, 1551, 1694, 2363, 2754, 2950)
#'t2 <- c(17, 42, 44, 48, 60, 72, 74, 95, 103, 108, 122, 144, 167, 170,
#'        183, 185, 193, 195, 197, 208, 234, 235, 254, 307, 315, 401, 445, 464, 484,
#'        528, 542, 547, 577, 580, 795, 855, 1366, 1577, 2060, 2412, 2486, 2796, 2802, 2934, 2988)
#'c1 <- c(rep(1, 43), 0, 0)
#'c2 <- c(rep(1, 39), rep(0, 6))
#'control_median_test(t1, c1, t2, c2, R = 500)
#'
#' @references 
#' Li, G., Tiwari, R.C., and Wells, M. (1996). "Quantile Comparison Functions in Two-Sample Problems: With Applications to Comparisons of Diagnostic Markers." Journal of the American Statistical Association, 91, 689-698.
#' 
#' Chakraborti, S., and Mukerjee, R. (1989), "A Confidence Interval for a Measure Associated With the Comparison of a Treatment With a Control," South African Statistical Journal, 23, 219-230.
#' 
#' Gastwirth, J. L., and Wang, J. L. (1988), "Control Percentile Test for Censored Data," Journal of Statistical Planning and Inference, 18, 267-276.
control_median_test <- function(t1, c1, t2, c2, R = 1000, seed = 1234){
  
  set.seed(seed)
  fit1 <- survfit(Surv(t1,c1)~1, conf.type = 'none')
  fit2 <- survfit(Surv(t2,c2)~1, conf.type = 'none')
  F.inv <- unname(quantile(fit1, prob = .5)) #Finv(p)
  G.inv <- unname(quantile(fit2, prob = .5)) #Ginv(p)
  
  Qp <- function(t1, c1, t2, c2){
    fit1 <- survfit(Surv(t1,c1)~1, conf.type = 'none') # Control
    fit2 <- survfit(Surv(t2,c2)~1, conf.type = 'none') # Experiment
    F.inv <- unname(quantile(fit1, prob = .5)) #Finv(p)
    G.inv <- unname(quantile(fit2, prob = .5)) #Ginv(p)
    
    if(is.na(F.inv)){
      warning("Median survival time could not be estimated for sample 1. Largest observed survival time was used as median.")
      F.inv <- max(t1)
    }
    
    if(is.na(G.inv)){
      warning("Median survival time could not be estimated for sample 2. Largest observed survival time was used as median.")
      G.inv <- max(t2)
    }
    
    G <- stepfun(fit2$time, c(1, fit2$surv))  #G
    Q <- G(F.inv)                             #Q = G[Finv(p)]
    return(Q)
  }
  
  Q <- Qp(t1, c1, t2, c2)
  
  #Bootstrapping:
  bootstrap.est <- numeric(R)
  for(i in 1:R){
    btsp1 <- sample(c(1:length(t1)), replace = TRUE) #bootstrapping for sample 1
    t1.bootstrap <- t1[btsp1]
    c1.bootstrap <- c1[btsp1]
    btsp2 <- sample(c(1:length(t2)), replace = TRUE) #bootstrapping for sample 2
    t2.bootstrap <- t2[btsp2]
    c2.bootstrap <- c2[btsp2]
    bootstrap.est[i] <- Qp(t1.bootstrap, c1.bootstrap, t2.bootstrap, c2.bootstrap)
  }
  se <- sd(bootstrap.est)
  Z = abs(Q-.5)/se
  pval <- 2*(1-pnorm(Z))
  
  plot(fit1, col = 'red', ylab = 'Estimated Survival Function', xlab = 'Time', main = 'Kaplan-Meier Estimates')
  lines(fit2, lty = 2, col = 'blue')
  legend('topright', c("KM-Estimate for Control Group", "KM-Estimate for Trt. Group"), lty = c(1, 2), col = c('red','blue'), bty = 'n', cex = .8)
  list("Med. for Control Group" = F.inv, "Med. for Trt. Group" = G.inv, "Z-score" = Z, "two-sided p-value" = pval)
}