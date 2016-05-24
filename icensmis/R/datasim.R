#' Simulate data including multiple outcomes from error-prone diagnostic tests 
#' or self-reports
#' 
#' @description This function simulates a data of N subjects with misclassified 
#'   outcomes, assuming each subject receives a sequence of pre-scheduled tests 
#'   for disease status ascertainment. Each test is subject to error, 
#'   characterized by sensitivity and specificity. An exponential distribution 
#'   is assumed for the time to event of interest. Three kinds of covariate 
#'   settings can be generated: one sample setting, two group setting, and 
#'   continuous covariates setting with each covariate sampled from i.i.d. N(0, 
#'   1). Two missing mechanisms can be assumed, namely MCAR and NTFP. The MCAR 
#'   setting assumes that each test is subject to a constant, independent 
#'   probability of missingness. The NTFP mechanism includes two types of 
#'   missingness - (1) incorporates a constant, independent, probability of 
#'   missing for each test prior to the first positive test result; and (2) all 
#'   test results after first positive are missing. The simulated data is in 
#'   longitudinal form with one row per test time.
#'   
#'   Covariate values, by default, are assumed to be constant. However, this 
#'   function can simulate a special case of time varying covariates. Under time
#'   varying covariates setting, each subject is assumed to have a change time 
#'   point, which is sampled from the visit times. We assume that each subject 
#'   has two sets of covariate values. Before his change time point, the 
#'   covariate values take from the first set, and second set after change time 
#'   point. Thus, each subject's distribution of survival time is two-piece 
#'   exponential distribution with different hazard rates.
#'   
#' @export
#' 
#' @param N total number of subjects to be simulated
#' @param blambda baseline hazard rate
#' @param testtimes a vector of pre-scheduled test times
#' @param sensitivity the sensitivity of test
#' @param specificity the specificity of test
#' @param betas a vector of regression coefficients of the same length as the 
#'   covariate vector. If betas = NULL then the simulated dataset corresponds to
#'   the one sample setting. If betas != NULL and twogroup != NULL then the 
#'   simulated dataset corresponds to the two group setting, and the first value
#'   of betas is used as the coefficient for the treatment group indicator. If 
#'   betas != NULL and twogroup = NULL, then the covariates are ~ i.i.d. N(0, 
#'   1), and the number of covariates is determined by the length of betas.
#' @param twogroup corresponds to the proportion of subjects allocated to the 
#'   baseline (reference) group in the two-group setting. For the two-group 
#'   setting, this variable should be between 0 and 1. For the one sample and 
#'   multiple (>= 2) covariate setting, this variable should be set to NULL. 
#'   That is, when betas !=NULL, set twogroup to equal the proportion of the 
#'   subjects in the baseline group to obtain a simulated dataset corresponding 
#'   to the two-group setting. Else, set twogroup=NULL to obtain either the one 
#'   sample setting (betas=NULL) or continuous covariates (betas !=NULL).
#' @param pmiss a value or a vector (must have same length as testtimes) of the 
#'   probabilities of each test being randomly missing at each test time. If 
#'   pmiss is a single value, then each test is assumed to have an identical 
#'   probability of missingness.
#' @param pcensor a value or a vector (must have same length as testtimes) of
#' the probability of censoring at each visit, assuming censoring process
#' is independent on other missing mechanisms.
#' @param design missing mechanism: "MCAR" or "NTFP"
#' @param negpred baseline negative predictive value, i.e. the probability of being 
#'   truely disease free for those who were tested (reported) as disease free at
#'   baseline. If baseline screening test is perfect, then negpred = 1.
#' @param time.varying indicator whether fitting a time varying covariate model 
#'   or not
#'   
#' @return simulated longitudinal form data frame
#'   
#' @details To simulate the one sample setting data, set betas to be NULL. To 
#'   simulate the two group setting data, set twogroup to equal the proportion 
#'   of the subjects in the baseline group and set betas to equal the 
#'   coefficient corresponding to the treatment group indicator(i.e. beta equals
#'   the log hazard ratio of the two groups). To simulate data with continuous 
#'   i.i.d. N(0, 1) covariates, set twogroup to be NULL and set betas to equal 
#'   the vector of coefficients of the covariates.
#'   
#' @examples
#' ## One sample setting  
#' simdata1 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = NULL, twogroup = NULL, pmiss = 0.3, design = "MCAR")
#' 
#' ## Two group setting, and the two groups have same sample sizes
#' simdata2 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = 0.7, twogroup = 0.5, pmiss = 0.3, design = "MCAR")
#'   
#' ## Three covariates with coefficients 0.5, 0.8, and 1.0
#' simdata3 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3,
#'   design = "MCAR", negpred = 1)
#' 
#' ## NTFP missing mechanism
#' simdata4 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3,
#'   design = "NTFP", negpred = 1)	 
#' 
#' ## Baseline misclassification
#' simdata5 <- datasim(N = 2000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3, 
#'   design = "MCAR", negpred = 0.97)  
#'   
#' ## Time varying covariates
#' simdata6 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3,
#'   design = "MCAR", negpred = 1, time.varying = TRUE)  
#'   
#' @useDynLib icensmis
#' @importFrom Rcpp evalCpp
#' @importFrom stats model.matrix optim pnorm qnorm rbinom rexp rnorm runif sd
#' 
datasim <- function(N, blambda, testtimes, sensitivity, specificity, 
                     betas = NULL, twogroup = NULL, pmiss = 0, pcensor = 0, design = "MCAR",
                     negpred = 1, time.varying = F) {
  
  #############################################################################
  # Input check
  ############################################################################# 
  ## type check
  stopifnot(is.numeric(N), is.numeric(N), is.numeric(testtimes),
            is.numeric(sensitivity), is.numeric(specificity), is.numeric(pmiss),
            is.numeric(negpred), is.logical(time.varying), is.character(design))
  ## range check
  stopifnot(length(N) == 1, N >= 1, is.finite(N), length(sensitivity) == 1, 
            sensitivity >= 0, sensitivity <= 1, length(specificity) == 1,
            specificity >= 0, specificity <= 1, length(negpred) == 1, negpred >= 0,
            negpred <= 1, all(testtimes > 0), all(is.finite(testtimes)),
            length(blambda) == 1, blambda > 0, is.finite(blambda),
            length(pmiss) == 1, pmiss < 1, pmiss >= 0, 
            design %in% c("MCAR", "NTFP"))
  if (!is.null(betas)) stopifnot(is.numeric(betas), is.finite(betas))
  if (!is.null(twogroup)) stopifnot(is.numeric(twogroup), length(twogroup) == 1, 
                                twogroup > 0, twogroup < 1)
  
  
  #############################################################################
  # Simulate covariates and event time
  #############################################################################
  testtimes <- sort(unique(testtimes))
  nbeta <- length(betas)
  ntest <- length(testtimes)
  ID <- rep(1:N, each = ntest)
  time <- rep(testtimes, times = N)
  if (!time.varying) {
    if (is.null(betas)) {
      covm <- NULL
      ET <- rexp(N, blambda)
    } else {
      if (!is.null(twogroup)) {
        nbeta <- 1
        covm <- matrix(rbinom(N, 1, 1 - twogroup), ncol = 1)
        colnames(covm) <- "group"
        lambda <- blambda*exp(c(covm*betas[1]))
      } else {
        covm <- matrix(rnorm(nbeta * N), ncol = nbeta)
        colnames(covm) <- paste0("cov", 1:nbeta)
        lambda <- blambda*exp(c(covm%*%betas))
      }
      ET <- rexp(N, lambda)     
      covm <- covm[ID, , drop = F]
    }     
  } else {
    if (is.null(betas)) stop("missing betas for time varying model")
    if (!is.null(twogroup)) 
      warning("two group setting is not supported for time varying model,
              twogroup = NULL assumed")
    chgtime <- sample(testtimes, N, replace = T)
    covm1 <- matrix(rnorm(nbeta*N), ncol = nbeta)
    covm2 <- matrix(rnorm(nbeta*N), ncol = nbeta)
    lambda1 <- blambda*exp(c(covm1%*%betas))
    lambda2 <- blambda*exp(c(covm2%*%betas))
    samplet <- function(ct, lam1, lam2) {
      x <- runif(1)
      cumhaz <- -log(1-x)
      ifelse(cumhaz < lam1*ct, cumhaz/lam1, (cumhaz - lam1*ct)/lam2 + ct)  
    }
    ET <- sapply(1:N, function(i) samplet(chgtime[i], lambda1[i], lambda2[i]))
    chgtime <- rep(chgtime, each = ntest)
    p1 <- time < chgtime
    p2 <- time >= chgtime
    covm <- p1*covm1[ID, , drop = F] + p2*covm2[ID, , drop = F]
    colnames(covm) <- colnames(covm1) <- paste0("cov", 1:nbeta)
  }
  ## baseline misclassification
  ET[rbinom(N, 1, 1 - negpred) == 1] <- 0
  ET <- ET[ID]
  
  #############################################################################
  # Simulate test result and apply missing
  #############################################################################  
  occur <- time > ET
  probs   <- ifelse(occur, sensitivity, 1 - specificity)
  result  <- rbinom(length(occur), 1, probs)
  if (is.null(covm)) {
    data <- data.frame(ID, testtime = time, result = result)
  } else {
    data <- data.frame(ID, covm, testtime = time, result = result)
  }
  notmiss <- rbinom(nrow(data), 1, 1 - pmiss)
  data <- data[notmiss == 1, ]
  
  #############################################################################
  # Apply NTFP
  ############################################################################# 
  if (design == "NTFP") {
    afterpos <- function(x) {
      npos <- cumsum(x == 1)
      (npos == 0) | (npos == 1 & x == 1)    
    }
    keep <- unlist(tapply(data$result, data$ID, afterpos))
    data <- data[keep, ]
  }
  
  #############################################################################
  # Apply censoring
  ############################################################################# 
  if (any(pcensor > 0)) {
    if (length(pcensor) == 1) pcensor <- rep(pcensor, ntest)
    pcensor <- c(pcensor, 1 - sum(pcensor))
    ctimes <- c(testtimes, Inf)
    data1 <- split(data, data$ID)
    data1 <- lapply(data1, function(x) {
      x$censortime <- sample(ctimes, 1, prob = pcensor)
      x
    })
    data1 <- unsplit(data1, data$ID)
    #data1 <- subset(data1, testtime < censortime)
    data1 <- data1[data1$testtime < data1$censortime, ]
    drop <- which(names(data1) == "censortime")
    data <- data1[, -drop]    
  }

  #############################################################################
  # Add baseline information for time-varying model
  ############################################################################# 
  if (time.varying) {
    base <- data.frame(ID = 1:N, covm1, testtime = 0, result = 0)
    base <- base[base$ID %in% unique(data$ID), ]
    data <- rbind(data, base)
    data <- data[order(data$ID, data$testtime), ]
  }
  
  row.names(data) <- NULL
  data
}