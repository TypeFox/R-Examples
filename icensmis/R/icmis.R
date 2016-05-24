#' Maximum likelihood estimation for settings of error-prone diagnostic tests 
#' and self-reported outcomes
#' 
#' This function estimates the baseline survival function evaluated at each test
#' time in the presence of error-prone diagnostic tests and self-reported 
#' outcomes. If there are covariates included in the dataset, it also estimates 
#' their coefficients assuming proportional hazards. The covariate values can be
#' either time independent or time varying The function can also be used to
#' incorporate misclassification of disease status at baseline (due to an
#' error-prone diagnostic procedure).
#' 
#' @param subject variable in data for subject id.
#' @param testtime variable in data for test time. Assume all test times are 
#'   non-negative. testtime = 0 refers to baseline visit (only used/needed if 
#'   the model is time varying covarites)
#' @param result variable in data for test result.
#' @param data the data to analyze.
#' @param sensitivity the sensitivity of test.
#' @param specificity the specificity of test.
#' @param formula a formula to specify what covariates to be included in the 
#'   model. If there is no covariate or one sample setting, set it to NULL. 
#'   Otherwise, input like ~x1 + x2 + factor(x3).
#' @param negpred baseline negative predictive value, i.e. the probability of 
#'   being truely disease free for those who were tested (reported) as disease 
#'   free at baseline. If baseline screening test is perfect, then negpred = 1.
#' @param time.varying indicator whether fitting a time varying covariate model 
#'   or not.
#' @param betai a vector of initial values for the regression coefficients 
#'   corresponding to the vector of covariates. If betai=NULL, then 0s are used 
#'   for the initial values. Otherwise, the length of betai must equal the 
#'   number of covariates.
#' @param initsurv initial value for survival function of baseline group in the 
#'   last visit time. It is used to compute initival values for survival 
#'   function at all visit times.
#' @param param parameterization for survival function used for optimization, 
#'   taking values 1, 2, or 3. There are 3 parameterizations available. param = 
#'   1: this parameterization uses the change in cumulative incidence in time 
#'   period j for baseline group as parameters, i.e. log(S[j]) - log(S[j+1]). 
#'   param = 2: simply use log of the parameters in param = 1 so that those 
#'   parameters are unbounded. param = 3: the first element is 
#'   log(-log(S[tau_1])) corresponding to log-log transformation of survival 
#'   function at first visit, while other parameters are corresponding to the 
#'   change in log-log of surival function, log(-log(S[j])) - log(-log(S[j-1])).
#'   In most cases, all parameters yield same results , while in some situations
#'   especially when two visit times are estimated to have same survival 
#'   functions, they may differ. Choose the one that works best (check 
#'   likelihood function)
#' @param ... other arguments passed to \code{\link{optim}} function. For 
#'   example, if the optimization does not converge, we can increase maxit in 
#'   the optim function.
#'   
#' @details The input data should be in longitudinal form with one row per test 
#'   time. Use \code{\link{datasim}} to simulate a dataset to see the sample 
#'   data structure. If time varying model is to be fitted, the baseline visit 
#'   must be provided so that the baseline covariate information can be 
#'   extracted. If an error is generated due to the optimization procedure, then
#'   we recommend trying different initial values.
#'   
#'   This likelihood-based approach is a function of the survival function
#'   evaluated at each unique test time in the dataset and the vector of
#'   regression coefficients as model parameters. Therefore, it works best for
#'   situations where there is a limited number of unique test times in the
#'   dataset. If there are a large number of unique test times, one solution is
#'   to group several test times together.
#'   
#' @export
#' 
#' @return A list of fitting results is returned with log-likelihood, estimated 
#'   coefficiets, estimated survival function, and estimated covariance matrix 
#'   for covariates.
#'   
#' @examples
#' ## One sample setting
#' simdata1 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = NULL, twogroup = NULL, pmiss = 0.3, design = "MCAR")
#' fit1 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata1,
#'   sensitivity = 0.7, specificity= 0.98, formula = NULL, negpred = 1)	
#'              			  
#' ## Two group setting, and the two groups have same sample sizes
#' simdata2 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = 0.7, twogroup = 0.5, pmiss = 0.3, design = "MCAR")					  
#' fit2 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata2,
#'   sensitivity = 0.7, specificity= 0.98, formula = ~group)						  
#'  
#' ## Three covariates with coefficients 0.5, 0.8, and 1.0
#' simdata3 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3,
#'   design = "MCAR", negpred = 1)					  
#' fit3 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata3,
#'   sensitivity = 0.7, specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 1)
#'  
#' ## Fit data with NTFP missing mechanism (the fitting is same as MCAR data)
#' simdata4 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3, 
#'   design = "NTFP", negpred = 1)						  
#' fit4 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata4,
#'   sensitivity = 0.7, specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 1)
#'               					  
#' ## Fit data with baseline misclassification
#' simdata5 <- datasim(N = 2000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3, 
#'   design = "MCAR", negpred = 0.97)						  
#' fit5 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata5,
#'   sensitivity = 0.7, specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 0.97)
#'   
#' ## Fit data with time varying covariates
#' simdata6 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7,
#'   specificity = 0.98, betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3,
#'   design = "MCAR", negpred = 1, time.varying = TRUE) 
#' fit6 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata6,
#'    sensitivity = 0.7, specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 1,
#'    time.varying = TRUE)  

icmis <- function(subject, testtime, result, data, sensitivity, specificity,
                   formula = NULL, negpred = 1, time.varying = F, 
                   betai = NULL, initsurv = 0.5, param = 1, ...){
  
  #############################################################################
  # Data Pre-processing: sort by id then by time
  #############################################################################
  id <- eval(substitute(subject), data, parent.frame())
  time <- eval(substitute(testtime), data, parent.frame())
  result <- eval(substitute(result), data, parent.frame())  
  ord <- order(id, time)
  if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    result <- result[ord]
    data <- data[ord, ]
  }  
  utime <- sort(unique(time))
  
  #############################################################################
  # Check variables and input parameters check
  ############################################################################# 
  stopifnot(is.numeric(sensitivity), is.numeric(specificity),
            is.numeric(negpred), is.logical(time.varying), is.numeric(initsurv))
  stopifnot(length(sensitivity) == 1, sensitivity >= 0, sensitivity <= 1, 
            length(specificity) == 1, specificity >= 0, specificity <= 1,
            length(negpred) == 1, negpred >= 0, negpred <= 1)
  stopifnot(length(initsurv) == 1, initsurv > 0, initsurv < 1)
  if (!all(result %in% c(0, 1))) stop("result must be 0 or 1")
  if (any(tapply(time, id, anyDuplicated)))
    stop("existing duplicated visit times for some subjects")  
  if (!all(utime >= 0 & utime < Inf))
    stop("existing negative or infinite visit times")
  stopifnot(length(param) == 1, param %in% c(1, 2, 3))
  if (param == 3 && time.varying) 
    stop("parameterization 3 is not available for time varying model")
  
  #############################################################################
  # Compute D matrix
  #############################################################################  
  timen0 <- (time != 0)
  Dm <- dmat(id[timen0], time[timen0], result[timen0], sensitivity,
             specificity, negpred)
  J <- ncol(Dm) - 1
  nsub <- nrow(Dm)
  
  #############################################################################
  # Create initial values based on initsurv and param
  #############################################################################
  if (param == 1) {
    lami <- rep(-log(initsurv)/J, J)
    tosurv <- function(x) exp(-cumsum(x))
    lowlam <- rep(0, J)
  } else if (param == 2) {
    lami <- log(rep(-log(initsurv)/J, J))
    tosurv <- function(x) exp(-cumsum(exp(x)))
  } else if (param == 3) {
    lami <- log(-log(seq(1, initsurv, length.out = J + 1)[-1]))
    lami <- c(lami[1], diff(lami))
    tosurv <- function(x) exp(-exp(cumsum(x)))
    lowlam <- c(-Inf, rep(0, J - 1))
  }
  
  #############################################################################
  # Function to create output
  #############################################################################
  output <- function(q) {
    loglik <- -q$value
    lam <- q$par[1:J]
    surv <- tosurv(lam)
    survival <- data.frame(time = utime[utime!=0], surv = surv)
    if (!is.null(formula)) {
      cov <- as.matrix(solve(q$hessian)[-(1:J), -(1:J)])
      rownames(cov) <- colnames(cov) <- beta.nm
      beta.fit <- q$par[-(1:J)]
      beta.sd <- sqrt(diag(cov))
      beta.z <- beta.fit/beta.sd
      p.value <- 2*(1-pnorm(abs(beta.z)))
      coef <- data.frame(coefficient = beta.fit, SE = beta.sd, z = beta.z,
                         p.value = p.value)   
    } else {
      cov <- NA
      coef <- NA
    }
    list(loglik = loglik, coefficient = coef, survival = survival, beta.cov = cov,
         nsub = nsub)
  }

  #############################################################################
  # No-covariate model (one sample case)
  #############################################################################
  if (is.null(formula)) {
    if (param == 1) {
      q <- optim(lami, loglikA0, gradlikA0, lower = lowlam, Dm = Dm, 
                 method = "L-BFGS-B", ...)
    } else if (param == 2) {
      q <- optim(lami, loglikB0, gradlikB0, Dm = Dm, method = "BFGS", ...)
    } else if (param == 3) {
      q <- optim(lami, loglikC0, gradlikC0, lower = lowlam, Dm = Dm, 
                 method = "L-BFGS-B", ...)
    }
    if (q$convergence != 0) 
      stop(paste("Not converged, code:", q$convergence)) 
    return(output(q))
  }
  
  #############################################################################
  # With-covariate model, time-independent covariates
  #############################################################################
  Xmat <- model.matrix(formula, data = data)[, -1, drop = F]
  if (nrow(Xmat) < nrow(data)) stop("missing values in used covariates")
  beta.nm <- colnames(Xmat)
  nbeta <- ncol(Xmat)
  if (is.null(betai)) {
    parmi <- c(lami, rep(0, nbeta))
  } else {
    if (length(betai) != nbeta) stop("length of betai does match number of beta")
    parmi <- c(lami, betai)
  }
  if (!time.varying) {
    uid <- getrids(id, nsub)
    Xmat <- Xmat[uid, , drop = F]
    if (param == 1) {
      q <- optim(parmi, loglikA, gradlikA, lower = c(lowlam, rep(-Inf, nbeta)),
                 Dm = Dm, Xmat = Xmat, method = "L-BFGS-B", hessian = T, ...)
    } else if (param == 2) {
      q <- optim(parmi, loglikB, gradlikB, Dm = Dm, Xmat = Xmat, method = "BFGS",
                 hessian = T, ...)
    } else if (param == 3) {
      q <- optim(parmi, loglikC, gradlikC, lower = c(lowlam, rep(-Inf, nbeta)),
                 Dm = Dm, Xmat = Xmat, method = "L-BFGS-B", hessian = T, ...)
    }    
    if (q$convergence != 0) 
      stop(paste("Not converged, code:", q$convergence)) 
    return(output(q))
  }
  
  #############################################################################
  # With-covariate model, time-varying covariates model
  #############################################################################  
  if (length(unique(id)) != length(unique(id[time==0]))) 
    stop("some subject(s) miss baseline information (time=0) 
         for time-varying model")
  TXmat <- timeMat(nsub, J+1, time, utime, Xmat)
  
  if (param == 1) {
    q <- optim(parmi, loglikTA, gradlikTA, lower = c(lowlam, rep(-Inf, nbeta)),
               Dm = Dm, TXmat = TXmat, method = "L-BFGS-B", hessian = T, ...)
  } else if (param == 2) {
    q <- optim(parmi, loglikTB, gradlikTB, Dm = Dm, TXmat = TXmat, method = "BFGS",
               hessian = T, ...)
  }  
  if (q$convergence != 0) 
    stop(paste("Not converged, code:", q$convergence)) 
  return(output(q))
}
