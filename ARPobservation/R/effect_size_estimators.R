# Generic error-checking function common to all of these functions

phase_validation <- function(observations, phase, base_level) {

  #recasts phase as a factor if it is just a vector
  phase <- factor(phase)
  
  if(!any(levels(phase) == base_level)){
    stop('Specified base level does not match a level in the data.')
  }
  
  #Checks to be sure there are only two levels
  if(length(levels(phase)) != 2){
    stop('This function requires exactly two levels in the phase variable')
  }
  
  if(length(phase) != length(observations)){
    stop('The length of \'observations\' and the length of \'phase\' need to be the same')
  }
  
  treat_level <- levels(phase)[(base_level != levels(phase))]
  
  return(c(base_level, treat_level))
}

#' @title Calculate log-response ratio, variance, and confidence interval
#' 
#' @description Estimates the log-response ratio (with or without bias correction),
#' the variance of the log-response ratio, and the confidence interval for a given
#' confidence level.
#' 
#' @param observations  Vector of observations
#' @param phase         Factor or vector indicating levels of the PIR measurements.
#' @param base_level a character string or value indicating the name of the baseline level.
#' @param conf_level Desired coverage rate of the calculated confidence interval. Default is \code{.95}.
#' @param bias_correct  Logical value indicating if the bias-corrected log-response ratio should be used. Default is \code{TRUE}
#' @param exponentiate  Logical value indicating if the log-respones ratio should be exponentiated.
#' 
#' @details The \code{observations} vector can be in any order corresponding to the factor or vector \code{phase}. 
#' The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". 
#' If there are more than two levels in \code{phase} this function will not work.
#' A value for \code{base_level} must be specified - if it is a chaaracter string it is case sensitive.
#' If \code{exponentiate = TRUE}, the log-ratio and the confidence interval will be exponentiated,
#' but the variance will be excluded from the output.
#' 
#' @return If \code{exponentiate = FALSE}, a list with three named entries. 
#' The first entry, \code{lRR}, is the estimated log-response ratio. 
#' The second entry, \code{V_lRR}, is the estimated variance of the log-response ratio.
#' The third entry, \code{CI}, is a vector containing the endpoints of a confidence
#' interval of \code{conf_level} coverage rate.
#' 
#' If \code{exponentiate = TRUE}, a list with two named entries.
#' The first entry, \code{RR}, is the estimated response ratio.  
#' The second entry, \code{CI}, is a vector containing the endpoints of a confidence
#' interval of \code{conf_level} coverage rate.
#'
#' @examples 
#' 
#' # Estimate the log response ratio and its variance for Carl from Moes dataset
#' data(Moes)
#' with(subset(Moes, Case == "Carl"),
#' logRespRatio(observations = outcome, phase = Phase, base_level = "No Choice"))
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
#' @export

logRespRatio <- function(observations, phase, base_level, conf_level = .95, 
                         bias_correct = TRUE, exponentiate = FALSE) {
  
  level_labels <- phase_validation(observations = observations, phase = phase, base_level = base_level)
    
  # calculate summary statistics for both samples, sort so that base level is first
  nObs <- table(phase)[level_labels]
  means <- tapply(observations, phase, mean)[level_labels]
  variances <- tapply(observations, phase, var)[level_labels]
  
  if (!all(means > 0)) stop('The mean of one or both phases is at the floor of 0.')
  
  if (bias_correct == TRUE) {
    BC <- log(means) - variances / (2 * nObs * means^2)
    lRR <- as.numeric(BC[2] - BC[1])
  } else {
    lRR <- log(means[2]) - log(means[1])
  }
  
  V_lRR <- sum(variances / (nObs * means^2))
  
  CI <- lRR + c(-1, 1) * qnorm(1-(1-conf_level)/2) * sqrt(V_lRR)
  
  if(exponentiate){
    return(list(RR = exp(lRR), CI = exp(CI)))
  }else{
  return(list(lRR = lRR, V_lRR = V_lRR, CI = CI))
  }
}

#' @title Prevalence bounds and confidence interval
#' 
#' @description Calculates a bound for the log of the prevalence ratio of two samples (referred to as baseline and treatment)
#' based on partial interval recording (PIR) data, assuming that the behavior follows an Alternating Renewal Process.
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param base_level a character string or value indicating the name of the baseline level.
#' @param mu_L the lower limit on the mean event duration
#' @param active_length length of the active observation interval
#' @param intervals the number of intervals in the sample of observations. Default is \code{NA}.
#' @param conf_level Coverage rate of the confidence interval. Default is \code{.95}.
#' @param exponentiate Logical value indicating if the log of the bounds and the confidence interval should be exponentiated. Default is \code{FALSE}.
#' 
#' @details The prevalence ratio estimate is based on the assumptions that 
#' 1) the underlying behavior stream follows an Alternating Renewal Process and
#' 2) the average event duration is greater than \code{mu_L}.
#' 
#' The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. 
#' The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". 
#' If there are more than two levels in \code{phase} this function will not work. 
#' A value for \code{base_level} must be specified - if it is a chaaracter string it is case sensitive.
#' 
#' For all of the following variables, the function assumes that if a vector of values is provided they are constant across all observations and simply uses the first value in that vector.
#' 
#' \code{mu_L} is the lower limit on the mean event durations. This is a single value assumed to hold for both samples of behavior
#' 
#' \code{active_length} This is the total active observation length. If the intervals are 15 seconds long but 5 seconds of each interval is reserved for recording purposes, \code{active_length= 10}. Times are often in seconds, but can be in any time unit.
#' 
#' \code{intervals} is the number of intervals in the observations. This is a single value and is assumed to be constant across both samples and all observations. This value is only relevant if the mean of one of the samples is at the floor or ceiling of 0 or 1. In that case it will be used to truncate the sample mean. If the sample mean is at the floor or ceiling and no value for \code{intervals} is provided, the function will stop.
#' 
#' @return A list with three named entries. 
#' The first entry, \code{estimate_bounds}, contains the lower and upper bound for the estimate of the prevalence ratio. 
#' The second entry, \code{estimate_SE}, contains the standard error of the estimate.
#' The third entry, \code{estimate_CI}, contains the lower and upper bounds for the confidence interval of the prevalence ratio. 
#'
#' @examples 
#' # Estimate bounds on the prevalence ratio for Carl from Moes dataset
#' data(Moes)
#' with(subset(Moes, Case == "Carl"),
#'  prevalence_bounds(PIR = outcome, phase = Phase, base_level = "No Choice", 
#'  mu_L = 10, active_length = active_length, intervals = intervals))
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
#' @export

prevalence_bounds <- function(PIR, phase, base_level, mu_L, active_length, intervals = NA, conf_level = 0.95, exponentiate = FALSE) {
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  level_labels <- phase_validation(observations = PIR, phase = phase, base_level = base_level)
  
  # calculate summary statistics for both samples, sort so that base level is first
  nObs <- table(phase)[level_labels]
  means <- tapply(PIR, phase, mean)[level_labels]
  variances <- tapply(PIR, phase, var)[level_labels]
  
  #change vectors to single values
  mu_L <- mu_L[1]
  active_length <- active_length[1]
  intervals <- intervals[1]
  
  if ((!all(means > 0) | !all(means < 1)) & is.na(intervals))  stop('One of the means is at the floor or ceiling of 0 or 1 and no value for intervals has been provided to perform truncation')
  
  means <- ifelse(means == 0, 1/(nObs * intervals),ifelse(means == 1, 1 - (1/(nobs * intervals)),means))
  
  
  #Natural log of the ratio of the two means
  R <- log(means[2]) - log(means[1])
  
  h <- log(mu_L + active_length) - log(mu_L)
  
  #calculate lower and upper bounds
  lower_bound <- as.numeric(R - h)
  upper_bound <- as.numeric(R + h)
  
  #variance of the log ratio
  variance_R <- as.numeric((variances[1]/(nObs[1] * (means[1]^2))) + (variances[2]/(nObs[2] * (means[2]^2))))
  
  #calculate the CI
  lower_CI <- as.numeric(R - (h + (qnorm(1-(1-conf_level)/2)*sqrt(variance_R))))
  upper_CI <- as.numeric(R + (h + (qnorm(1-(1-conf_level)/2)*sqrt(variance_R))))
  
  #exponentiates the values, if desired
  if(exponentiate == TRUE){
    lower_bound <- exp(lower_bound)
    upper_bound <- exp(upper_bound)
    lower_CI <- exp(lower_CI)
    upper_CI <- exp(upper_CI)
  }
  
  return(list(estimate_bounds = c(lower_bound = lower_bound, upper_bound = upper_bound), 
              estimate_SE = sqrt(variance_R),
              estimate_CI = c(lower_CI = lower_CI, upper_CI = upper_CI)))
}

#' Incidence bounds and confidence interval
#' 
#' @description Calculates a bound for the log of the incidence ratio of two samples 
#' (referred to as baseline and treatment) based on partial interval recording (PIR) data, 
#' assuming that the behavior follows an Alternating Renewal Process.
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param base_level a character string or value indicating the name of the baseline level.
#' @param mu_U the upper limit on the mean event duration
#' @param p upper limit on the probability that the interim time between behavioral events is less than the active interval
#' @param active_length length of the active observation interval
#' @param intervals the number of intervals in the sample of observations. Default is \code{NA}
#' @param conf_level Coverage rate of the confidence interval. Default is \code{.95}.
#' @param exponentiate Logical value indicating if the log of the bounds and the confidence interval should be exponentiated. Default is \code{FALSE}.
#' 
#' @details The incidence ratio estimate is based on the assumptions that 
#' 1) the underlying behavior stream follows an Alternating Renewal Process,
#' 2) the average event duration is less than \code{mu_U}, and
#' 3) the probability of observing an interim time less than the active interval length is less than \code{p}.
#' 
#' The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", or "0" and "1". If there are more than two levels in \code{phase} this function will not work. A value for \code{base_level} must be specified - if it is a chaaracter string it is case sensitive.
#'  
#' For all of the following variables, the function assumes that if a vector of values is provided they are constant across all observations and simply uses the first value in that vector.
#' 
#' \code{mu_U} is the upper limit on the mean event durations. This is a single value assumed to hold for both samples of behavior
#' 
#' \code{active_length} This is the total active observation length. If the intervals are 15 seconds long but 5 seconds of each interval is reserved for recording purposes, \code{active_length= 10}. Times are often in seconds, but can be in any time unit.
#' 
#' \code{intervals} is the number of intervals in the observations. This is a single value and is assumed to be constant across both samples and all observations. This value is only relevant if the mean of one of the samples is at the floor or ceiling of 0 or 1. In that case it will be used to truncate the sample mean. If the sample mean is at the floor or ceiling and no value for \code{intervals} is provided, the function will stop.
#' 
#' @return A list containing two named vectors and a single named number. 
#' The first entry, \code{estimate_bounds}, contains the lower and upper bound for the estimate of the incidence ratio. 
#' The second entry, \code{estimate_SE}, contains the standard error of the estimate.
#' The third entry, \code{estimate_CI}, contains the lower and upper bounds for the confidence interval of the incidence ratio. 
#' 
#' @examples 
#' 
#' # Estimate bounds on the incidence ratio for Ahmad from the Dunlap dataset
#' data(Dunlap)
#' with(subset(Dunlap, Case == "Ahmad"),
#' incidence_bounds(PIR = outcome, phase = Phase, base_level = "No Choice", 
#'                  mu_U = 10, p = .15, active_length = active_length, intervals = intervals))
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
#' @export

incidence_bounds <- function(PIR, phase, base_level, mu_U, p, active_length, 
                             intervals = NA, conf_level = 0.95, 
                             exponentiate = FALSE) {
  
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  level_labels <- phase_validation(observations = PIR, phase = phase, base_level = base_level)
  
  # calculate summary statistics for both samples, sort so that base level is first
  nObs <- table(phase)[level_labels]
  means <- tapply(PIR, phase, mean)[level_labels]
  variances <- tapply(PIR, phase, var)[level_labels]
  
  mu_U <- mu_U[1]
  active_length <- active_length[1]
  intervals <- intervals[1]
  
  if ((!all(means > 0) | !all(means < 1)) & is.na(intervals))  stop('One of the means is at the floor or ceiling of 0 or 1 and no value for intervals has been provided to perform truncation')
  
  means <- ifelse(means == 0, 1/(nObs * intervals),ifelse(means == 1, 1 - (1/(nobs * intervals)),means))
  
  
  R <- log(means[2]) - log(means[1])
  
  h <- log(mu_U + active_length) - log(1-p) - log(active_length)
  
  lower_bound <- as.numeric(R - h)
  upper_bound <- as.numeric(R + h)
  
  variance_R <- as.numeric((variances[1]/(nObs[1] * (means[1]^2))) + (variances[2]/(nObs[2] * (means[2]^2))))
  
  lower_CI <- as.numeric(R - (h + (qnorm(1-(1-conf_level)/2) * sqrt(variance_R))))
  upper_CI <- as.numeric(R + (h + (qnorm(1-(1-conf_level)/2) * sqrt(variance_R))))
 
  if(exponentiate == TRUE){
    exp(lower_bound)
    exp(upper_bound)
    exp(lower_CI)
    exp(upper_CI)
  }
                   
  return(list(estimate_bounds = c(lower_bound = lower_bound,upper_bound = upper_bound), 
              estimate_SE = sqrt(variance_R),
              estimate_CI = c(lower_CI = lower_CI,upper_CI = upper_CI)))
}

#' @title Interim bounds and confidence interval
#' 
#' @description Calculates a bound for the log of the ratio of interim time of two samples 
#' (referred to as baseline and treatment) based on partial interval recording (PIR) data,
#' assuming that the average event durations are equal across samples and that 
#' interim times are exponentially distributed.
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param base_level a character string or value indicating the name of the baseline level.
#' @param conf_level Desired coverage rate of the calculated confidence interval. Default is \code{.95}.

#' @param intervals the number of intervals in the sample of observations. Default is \code{NA}
#' @param exponentiate Logical value indicating if the log of the bounds and the confidence interval should be exponentiated. Default is \code{FALSE}.
#' 
#' @details The interim ratio estimate is based on the assumptions that 
#' 1) the underlying behavior stream follows an Alternating Renewal Process,
#' 2) the average event durations in each sample are equal, and 
#' 3) interim times follow exponential distributions.
#' 
#' The \code{PIR} vector can be in any order corresponding to the factor or vector \code{phase}. 
#' The levels of \code{phase} can be any two levels, such as "A" and "B", "base" and "treat", 
#' or "0" and "1". If there are more than two levels in \code{phase} this function will not work. 
#' A value for \code{base_level} must be specified; if it is a chaaracter string it is case sensitive.
#' 
#' \code{intervals} is the number of intervals in the observations. 
#' This is a single value and is assumed to be constant across both samples and all observations. 
#' If intervals is sent as a vector instead of a single value, the first value in the vector will 
#' be used. This value is only relevant if the mean of one of the samples is at the floor or ceiling 
#' of 0 or 1. In that case it will be used to truncate the sample mean. If the sample mean is at the 
#' floor or ceiling and no value for \code{intervals} is provided, the function will stop.
#' 
#' @return A list with three named entries 
#' The first entry, \code{estimate_bounds}, contains the lower and upper bound for the estimate of the prevalence ratio. 
#' The second entry, \code{estimate_SE}, contains the standard errors for the upper and lower bounds. 
#' The third entry, \code{estimate_CI}, contains the lower and upper bounds for the confidence interval of the prevalence ratio. 
#' 
#' @examples 
#' # Estimate bounds on the interim time ratio for Carl from the Moes dataset
#' data(Moes)
#' with(subset(Moes, Case == "Carl"),
#' interim_bounds(PIR = outcome, phase = Phase, base_level = "No Choice"))
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
#' @export

interim_bounds <- function(PIR, phase, base_level, 
                           conf_level = 0.95, intervals = NA, exponentiate = FALSE) {
  
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  level_labels <- phase_validation(observations = PIR, phase = phase, base_level = base_level)
  
  # calculate summary statistics for both samples, sort so that base level is first
  nObs <- table(phase)[level_labels]
  means <- tapply(PIR, phase, mean)[level_labels]
  variances <- tapply(PIR, phase, var)[level_labels]
  
  intervals <- intervals[1]
  
  if ((!all(means > 0) | !all(means < 1)) & is.na(intervals))  stop('One of the means is at the floor or ceiling of 0 or 1 and no value for intervals has been provided to perform truncation')
  
  means <- ifelse(means == 0, 1/(nObs * intervals),ifelse(means == 1, 1 - (1/(nobs * intervals)),means))
  
  #the logit and complimentary log-log functions
  logit <- function(x) {log(x) - log(1-x)}
  cll <- function(x) {log(-1 * log(1-x))}
  
  #This code checks to see which transformation is the appropriate value for
  #the lower and upper bounds
  if(means[1] > means[2]){
    f_lower <- as.numeric(logit(means[2]) - logit(means[1]))
    f_upper <- as.numeric(cll(means[2]) - cll(means[1]))
  }else{
    f_lower <- as.numeric(cll(means[2]) - cll(means[1]))
    f_upper <- as.numeric(logit(means[2]) - logit(means[1]))
  }
  
  
  #This is the variance of the log ratio of the sample means
  var_LOR <- as.numeric((variances[1]/(nObs[1] * (means[1]^2) * (1- means[1])^2)) + 
             (variances[2]/(nObs[2] * (means[2]^2) * (1- means[2])^2)))
  
  #This is the variance of the complimentary log-log ratio
  var_CLR <- as.numeric((variances[1]/(nObs[1] * (1- means[1])^2 * (log(1-means[1]))^2)) + 
             (variances[2]/(nObs[2] * (1- means[2])^2 *(log(1-means[2]))^2)))
  
  #The ratio and z_conf is used in determining which variance is appropriate
  cll_ratio <- cll(means[2]) - cll(means[1])
  
  z_conf <- qnorm(1-(1-conf_level)/2)
  
  if(cll_ratio <= (z_conf * sqrt(var_LOR))){
    var_f_lower <- var_LOR
  } else{
    var_f_lower <- var_CLR
  }
  
  if(cll_ratio < -(z_conf * sqrt(var_LOR))){
    var_f_upper <-  var_CLR
  }else{
    var_f_upper <- var_LOR
  }
  
  lower_CI <- f_lower - (z_conf * sqrt(var_f_lower))
  upper_CI <- f_upper + (z_conf * sqrt(var_f_upper))
  
  if(exponentiate == TRUE)
  {
    f_lower <- exp(f_lower)
    f_upper <- exp(f_upper)
    lower_CI <- exp(lower_CI)
    upper_CI <- exp(upper_CI)
  }
  
  return(list(estimate_bounds = c(lower_bound = f_lower,upper_bound = f_upper),
              estimate_SE = c(lower_SE = sqrt(var_f_lower), 
                              upper_SE = sqrt(var_f_upper)),
              estimate_CI = c(lower_CI = lower_CI, upper_CI = upper_CI)))
}


# Zeta as a function of phi and E(Y)

PIR_Zeta <- function(phi, ExY, active) -1 * (1-phi) * log((ExY-1)/(phi-1)) / active

# Var(Y) as a function of phi and E(Y)

PIR_VarEx <- function(phi, ExY, active, L, K){
  
  zeta = PIR_Zeta(phi, ExY, active)
  
  if (phi == 0) {
    vary = ExY * (1- ExY) / K 
  } else {
    k = 1:(K-1) #vector for summation term
    vary = ((1/K) * ExY * (1-ExY)) * ( 1 + (2 * phi / K / ExY) * sum((K - k) * exp((zeta * active / phi)  - zeta * k * L /phi / (1- phi) / K)) )
  }  
  
  return(vary)
}

# Finds the root for phi using the variance and mean

PIR_Phi <- function(ExY, VarY, nObs, active, L, K) {
  fun <- function(phi) PIR_VarEx(phi, ExY, active, L, K) / VarY - 1
  uniroot(fun, interval = c(0, ExY), tol = .Machine$double.eps^0.5)$root
}

# DEPRECATED
# generates multiple samples of PIR data using the r_behavior_stream and interval_recording
# functions from ARPobservation. Useful for simulating or performing bootstraps

generatePIRData <- function(nObs, phi, zeta, active, K, rest = 0, iterations = 1) {
  
  #necessary to recast as numbers because PIR_MOM sends them as lists with one element
  phi <- as.numeric(phi) 
  zeta <- as.numeric(zeta)
  #calculating the means of the distributions of event durations and interim
  #times from estimates of prevalence and incidence
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  
  #the full interval length is the sum of the active and rest (recording) lengths
  intervalLength = active + rest
  
  #this block generates a matrix of PIR data with each column corresponding to a
  #a sample of n = nObs observations
  sampleObs <- replicate(iterations, {
    BS <- r_behavior_stream(n = nObs, mu = mu, lambda = lambda, 
                            F_event = F_exp(), F_interim = F_exp(), 
                            stream_length = intervalLength * K)
    interval_recording(BS = BS, interval_length = intervalLength, rest_length = rest)
  })
  
  mean <- colMeans(sampleObs)
  
  variance <- apply(sampleObs, 2, var)
  
  sampleData <- data.frame(mean = mean, variance = variance)
  
  return(sampleData)
}

# A packaging function that invokes the previous set of functions to neatly
# provide an estimate of phi and zeta for a sample or vector of samples

PIR_inv <- function(ExY, VarY, nObs, active, L, K, varTrunc = 1 / (nObs * K^2)) {
  
  ExY <- pmin(pmax(ExY, 1 / (nObs * K)), 1 - 1 / (nObs * K))
  VarY <- pmin(pmax(VarY, ExY * (1 - ExY) / K + varTrunc), ExY * (1 - ExY) - varTrunc)
  
  phi <- mapply(PIR_Phi, ExY = ExY, VarY = VarY, 
                nObs = nObs, active = active, L = L, K = K)
  
  zeta <- mapply(PIR_Zeta, phi = phi, ExY = ExY, active = active)
  
  ests <- data.frame(phi = phi, zeta = zeta)
  
  return (ests)
  
}


# bootstraps the confidence intervals for a pair of samples with estimates of phi and zeta

PIRbootstrappair <- function(nObs, phi, zeta, active, rest, K, iterations, alpha,
                            exponentiate, seed = NULL){
  # set seed if one is supplied
  if (!is.null(seed)) set.seed(seed)

  #the total length of the observation is the total interval length times the number of intervals
  L = (active + rest) * K
  
  #calculating the means of the distributions of event durations and interim
  #times from estimates of prevalence and incidence
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  
  # Simulate first set of parms
  sampleData0 <- replicate(n = iterations, r_PIR(n = nObs[1], mu = mu[1],
                                                 lambda = lambda[1],
                                                 stream_length = L,
                                                 F_event = F_exp(), 
                                                 F_interim = F_exp(),
                                                 interval_length = active + rest,
                                                 rest_length = rest, 
                                                 summarize = TRUE))
  
  mean0 <- colMeans(sampleData0)
  variance0 <- apply(sampleData0, 2, var)
  
  # Returns a dataframe of estimates of phi and zeta
  ests0 <- PIR_inv(ExY = mean0, VarY = variance0, nObs = nObs[1], 
                   active = active, K = K, L = L)
  
  # get the confidence interval bounds for the first sample
  pbounds0 <- quantile(ests0$phi, probs = c(alpha/2, 1-alpha/2))
  zbounds0 <- quantile(ests0$zeta, probs = c(alpha/2, 1-alpha/2))
  
  # simulate second set of parms
  sampleData1 <- replicate(n = iterations, r_PIR(n = nObs[2], mu = mu[2],
                                                 lambda = lambda[2],
                                                 stream_length = L,
                                                 F_event = F_exp(), 
                                                 F_interim = F_exp(),
                                                 interval_length = active + rest,
                                                 rest_length = rest, 
                                                 summarize = TRUE))
  
  mean1 <- colMeans(sampleData1)
  variance1 <- apply(sampleData1, 2, var)
  
  ests1 <- PIR_inv(ExY = mean1, VarY = variance1, nObs = nObs[2], 
                   active = active, K = K, L = L)
  
  pbounds1 <- quantile(ests1$phi, probs = c(alpha/2, 1-alpha/2))
  zbounds1 <- quantile(ests1$zeta, probs = c(alpha/2, 1-alpha/2))
  
  # calculate log ratio value and confidence interval bounds
  plogratio <- log(ests1$phi/ests0$phi)
  plogbounds <- quantile(plogratio, probs = c(alpha/2, 1-alpha/2))
  
  zlogratio <- log(ests1$zeta/ests0$zeta)
  zlogbounds <- quantile(zlogratio, probs = c(alpha/2, 1-alpha/2))
  
  if(exponentiate){
    pratbounds <- exp(plogbounds)
    zratbounds <- exp(zlogbounds)
    
    results <- cbind(phi = c(phi[1], phi[2], as.numeric(phi[2])/as.numeric(phi[1])),
                     phi_lower_CI = c(pbounds0[1], pbounds1[1], pratbounds[1]),
                     phi_upper_CI = c(pbounds0[2], pbounds1[2], pratbounds[2]),
                     zeta = c(zeta[1], zeta[2], as.numeric(zeta[2])/as.numeric(zeta[1])),
                     zeta_lower_CI = c(zbounds0[1], zbounds1[1], zratbounds[1]),
                     zeta_upper_CI = c(zbounds0[2], zbounds1[2], zratbounds[2]))
  }else{
    results <- cbind(phi = c(phi[1], phi[2], log(as.numeric(phi[2])/as.numeric(phi[1]))),
                     phi_lower_CI = c(pbounds0[1], pbounds1[1], plogbounds[1]),
                     phi_upper_CI = c(pbounds0[2], pbounds1[2], plogbounds[2]),
                     zeta = c(zeta[1], zeta[2], log(as.numeric(zeta[2])/as.numeric(zeta[1]))),
                     zeta_lower_CI = c(zbounds0[1], zbounds1[1], zlogbounds[1]),
                     zeta_upper_CI = c(zbounds0[2], zbounds1[2], zlogbounds[2]))
  }
  
  return(results)
  
}

#deprecated version of PIRbootstrappair that didn't use our general purpose
#data generating code

PIRbootstrappair_old <- function(nObs, phi, zeta, active, rest, K, iterations, alpha,
                             seed = NULL){
  # set seed if one is supplied
  if (!is.null(seed)) set.seed(seed)
  
  #the total length of the observation is the total interval length times the number of intervals
  L = (active + rest) * K
  
  # Simulate first set of parms
  sampleData0 <- generatePIRData(nObs = nObs[1], phi = phi[1], zeta = zeta[1],
                                 active = active, K = K, rest = rest, 
                                 iterations = iterations)
  
  
  # Returns a dataframe of estimates of phi and zeta
  ests0 <- with(sampleData0, PIR_inv(ExY = mean, VarY = variance, 
                                     nObs = nObs[1], active = active, K = K, 
                                     L = L))
  
  # get the confidence interval bounds for the first sample
  pbounds0 <- quantile(ests0$phi, probs = c(alpha/2, 1-alpha/2))
  zbounds0 <- quantile(ests0$zeta, probs = c(alpha/2, 1-alpha/2))
  
  # simulate second set of parms
  sampleData1 <- generatePIRData(nObs = nObs[2], phi = phi[2], zeta = zeta[2],
                                 active = active, K = K, rest = rest, 
                                 iterations = iterations)
  
  
  ests1 <- with(sampleData1, PIR_inv(ExY = mean, VarY = variance, 
                                     nObs = nObs[2], active = active, K = K, 
                                     L = L))
  pbounds1 <- quantile(ests1$phi, probs = c(alpha/2, 1-alpha/2))
  zbounds1 <- quantile(ests1$zeta, probs = c(alpha/2, 1-alpha/2))
  
  # calculate log ratio value and confidence interval bounds
  plogratio <- log(ests1$phi/ests0$phi)
  plogbounds <- quantile(plogratio, probs = c(alpha/2, 1-alpha/2))
  
  zlogratio <- log(ests1$zeta/ests0$zeta)
  zlogbounds <- quantile(zlogratio, probs = c(alpha/2, 1-alpha/2))
  
  results <- cbind(phi = c(phi[1], phi[2], log(as.numeric(phi[2])/as.numeric(phi[1]))),
                   phi_lower_CI = c(pbounds0[1], pbounds1[1], plogbounds[1]),
                   phi_upper_CI = c(pbounds0[2], pbounds1[2], plogbounds[2]),
                   zeta = c(zeta[1], zeta[2], log(as.numeric(zeta[2])/as.numeric(zeta[1]))),
                   zeta_lower_CI = c(zbounds0[1], zbounds1[1], zlogbounds[1]),
                   zeta_upper_CI = c(zbounds0[2], zbounds1[2], zlogbounds[2]))
  
  return(results)
  
}

#' @title Moment estimator for prevalence and incidence, with bootstrap confidence intervals
#' 
#' @description Estimates prevalance and incidence for two samples, 
#' along with the ratios of each parameter, assuming that the behavior follows 
#' an `. Also provides boostrap confidence intervals.
#' 
#' @param PIR vector of PIR measurements
#' @param phase factor or vector indicating levels of the PIR measurements.
#' @param base_level a character string or value indicating the name of the baseline level.
#' @param intervals the number of intervals in the sample of observations
#' @param interval_length the total length of each interval
#' @param rest_length length of the portion of the interval devoted to recording. Default is \code{0}
#' @param Bootstraps desired number of bootstrap replicates. Default is \code{2000}
#' @param conf_level Desired coverage rate of the calculated confidence interval. Default is \code{.95}.
#' @param exponentiate a logical indicating whether the row corresponding to the ratio of treatment to baseline should be exponentiated, with the default as \code{FALSE}.
#' @param seed seed value set in order to make bootstrap results reproducible. Default is \code{null}
#' 
#' @details The moment estimators are based on the assumption that the 
#' underlying behavior stream follows an Alternating Poisson Process, in which both the
#' event durations and interim times are exponentially distributed.
#' 
#' @details The \code{PIR} vector can be in any order corresponding to the factor or vector 
#' \code{phase}. The levels of \code{phase} can be any two levels, such as "A" and "B", 
#' "base" and "treat", or "0" and "1". If there are more than two levels in \code{phase}
#' this function will not work. A value for \code{base_level} must be specified - if it is a 
#' character string it is case sensitive.
#' 
#' \code{intervals}, \code{interval_length}, and \code{rest_length} are all single values that 
#' are assumed to be held constant across both samples and all observation sessions. 
#' If vectors of values are provided for these variables, it is assumed that the first value 
#' in each vector is constant across all observations.
#' 
#' \code{interval_length} This is the total length of each individual interval. 
#' Sometimes a portion of the interval is set aside for recording purposes, in which case 
#' \code{rest_length} should be set to the length of time devoted to recording. 
#' The default assumption is that there is no recording time. 
#' The length of time devoted to active observation is calculated as \code{interval_length - rest_length}.
#' 
#' At the default setting of \code{bootstraps = 2000}, PIR_MOM takes just under six seconds to run on an Intel Core i5-2410M processor.
#' 
#' @return A dataframe with six columns and three rows corresponding to baseline, treatment, 
#' and the log ratio or ratio (depending upon the value of \code{exponentiate}) of treatment to baseline
#' 
#' @examples 
#' 
#' # Estimate prevalence and incidence ratios for Carl from the Moes dataset
#' data(Moes)
#' with(subset(Moes, Case == "Carl"),
#' PIR_MOM(PIR = outcome, phase = Phase, intervals = intervals, 
#' interval_length = (active_length + rest_length), rest_length = rest_length, 
#' base_level = "No Choice", seed = 149568373))
#' 
#' @author Daniel Swan <dswan@@utexas.edu>
#' @export

PIR_MOM <- function(PIR, phase, base_level, intervals, interval_length, rest_length = 0, Bootstraps = 2000, conf_level = 0.95, exponentiate = FALSE, seed = NULL) {
  
  if(length(which(PIR > 1 | PIR < 0)) > 0 | sum(is.na(PIR)) > 0) {
    stop('Values for PIR must be between 0 and 1 and cannot be NA')
  }
  
  level_labels <- phase_validation(observations = PIR, phase = phase, base_level = base_level)
  
  # calculate summary statistics for both samples, sort so that base level is first
  nObs <- table(phase)[level_labels]
  means <- tapply(PIR, phase, mean)[level_labels]
  variances <- tapply(PIR, phase, var)[level_labels]
  
  intervals <- intervals[1]
  interval_length <- interval_length[1]
  rest_length <- rest_length[1]
  
  if ((!all(means > 0) | !all(means < 1)) & is.na(intervals))  stop('One of the means is at the floor or ceiling of 0 or 1 and no value for intervals has been provided to perform truncation')
  
  #Get an estimate of phi and zeta using the moment estimator for baseline and treatment phases
  base_ests <- PIR_inv(ExY = means[1], VarY = variances[1],
                       nObs = nObs[1], active = interval_length - rest_length,
                       L = intervals * interval_length, K = intervals)  
  
  treat_ests <- PIR_inv(ExY = means[2], VarY = variances[2],
                      nObs = nObs[2], active = interval_length - rest_length, 
                      L = intervals * interval_length, K = intervals)
  
  
  #Bootstrap the confidence intervals
  results <- PIRbootstrappair(nObs = c(nObs[[1]], nObs[[2]]),
                   phi = c(base_ests[[1]], treat_ests[[1]]),
                   zeta = c(base_ests[[2]], treat_ests[[2]]),
                   active = interval_length - rest_length,
                   rest = rest_length,
                   K = intervals,
                   iterations = Bootstraps,
                   alpha = 1-conf_level,
                   exponentiate = exponentiate,
                   seed = seed)
  
  #Set the row names appropriately based on the levels in the "phase" variable
  if(exponentiate){
    row.names(results) <- c(base_level, 
                            levels(phase)[which(levels(phase) != base_level)], 
                            paste(levels(phase)[which(levels(phase) != base_level)], 
                            base_level, sep = "/"))
  }else{
    row.names(results) <- c(base_level, 
                            levels(phase)[which(levels(phase) != base_level)], 
                            paste0("log(",paste(levels(phase)[which(levels(phase) != base_level)], 
                            base_level, sep = "/"),")"))
  }
  return(results)
}

#' Dunlap et al.(1994) data
#' 
#' Single case design data measured with partial interval recording from a study of the effect of providing Choice between academic activities on the disruptive behavior of three elementary school students with emotional and behavioral disorders. For this data "No Choice" is the baseline phase. Data were extracted from the figures in the publicaton.
#' 
#' @usage Dunlap
#' 
#' @format A data frame with 58 observations on 7 variables
#' 
#' \itemize{
#' \item [,1] \code{Case} The participant for whom the observation took place
#' \item [,2] \code{Phase} The level of the observation ("Choice" vs. "No Choice")
#' \item [,3] \code{Session} The observation session # for each participant
#' \item [,4] \code{outcome} The summary PIR measurement for the observation session
#' \item [,5] \code{active_length} The length of the active observation interval, in seconds
#' \item [,6] \code{rest_length} The length of the recording interval, in seconds
#' \item [,7] \code{intervals} The total number of intervals in the observation session 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 53940 rows and 10 variables
#' @name Dunlap
#' @references 
#' Dunlap, G., DePerczel, M., Clarke, S., Wilson, D., Wright, S., White, R., & Gomez, A. (1994). Choice making to promote adaptive behavior for students with emotional and behavioral challenges. Journal of Applied Behavior Analysis, 27 (3), 505-518.
NULL

#' Moes(1998) data
#' 
#' @description Single-case design data from a study that using partial interval recording (PIR) examining the impact of choice-making in a homework tutoring context on disruptive behavior. In this data "No Choice" is the baseline phase. Data were extracted from the figure in the publication.
#' @usage Moes
#' 
#' @format A data frame with 80 observations on 7 variables
#' 
#' \itemize{
#' \item [,1] \code{Case} The participant for whom the observation took place
#' \item [,2] \code{Phase} The level of the observation ("Choice" vs. "No Choice")
#' \item [,3] \code{Session} The observation session # for each participant
#' \item [,4] \code{outcome} The summary PIR measurement for the observation session
#' \item [,5] \code{active_length} The length of the active observation interval, in seconds
#' \item [,6] \code{rest_length} The length of the recording interval, in seconds
#' \item [,7] \code{intervals} The total number of intervals in the observation session 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 53940 rows and 10 variables
#' @name Moes
#' @references 
#' Moes, D. R. (1998). Integrating choice-making opportunities within teacher-assigned academic tasks to facilitate the performance of children with autism. Research and Practice for Persons with Severe Disabilities, 23 (4), 319-328.
NULL
