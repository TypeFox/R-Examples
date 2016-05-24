#----------------------------------------------------------------------#
# Test and identify a subgroup with an enhanced treatment effect       #
#----------------------------------------------------------------------#
#                                                                      #
#  Inputs                                                              #
#                                                                      #
# outcome  A formula object.                                           #
#          The linear model for the outcome regression.                #
#          The left-hand-side variable must be the response.           #
#          R function lm will be used to estimate model parameters.    #
#                                                                      #
# propen   A formula object.                                           #
#          The model for the propensity for treatment.                 #
#          The left-hand-side variable must be the treatment variable. #
#          R function glm will be used with input option               #
#          family = binomial(link="logit") to estimate model parameters#
#                                                                      #
# data     A data.frame object.                                        #
#          All covariates, treatment, and response variables.          #
#          Note that the treatment must be binary and that the         #
#          response must be continuous.                                #
#                                                                      #
# K        An integer object.                                          #
#          The number of random sampled points on the unit ball        #
#          surface \eqn{\{\theta:||\theta||^2=1\}}.                    #
#          These randomly sampled points are used for approximating    #
#          the Gaussian process in null and local alternative          #
#          distributions of the test statistic with multivariate       #
#          normal distributions.                                       #
#          It is recommended that K be set to 10^p, where p is the     #
#          number of parameters in theta0.                             #
#          Note that it is recommended that the number of covariates   #
#          be less than 10 for this implementation.                    #
#          Default value is 1000.                                      #
#                                                                      #
# M        An integer object.                                          #
#          The number of resamplings of the perturbed test statistic.  #
#          This sample is used to calculate the critical value of the  #
#          test.                                                       #
#          Default and minimum values are 1000.                        #
#                                                                      #
# seed     An integer object or NULL.                                  #
#          If set, the seed for generating random points on the unit   #
#          ball surface.                                               #
#          If NULL, current seed in R environment is used.             #
#                                                                      #
#  Outputs                                                             #
#                                                                      #
#  A list is returned containing.                                      #
#                                                                      #
#    outcome  An lm object.                                            #
#             The object returned by the lm fit of the outcome.        #
#                                                                      #
#    propen   A glm object.                                            #
#             The object returned by the glm fit of the propensity.    #
#                                                                      #
#    p_value  A numeric object.                                        #
#             The p-value of the test.                                 #
#                                                                      #
#    theta    A named numeric vector.                                  #
#             The change-plane parameter estimates for subgroup.       #
#                                                                      #
#    prop     A numeric object.                                        #
#             The proportion of sampled points on \eqn{\theta} unit    #
#             ball surface that are used for calculating test          #
#             statistic. When calculating the test statistic, the      #
#             subgroup indicated by some values of \eqn{\theta}        #
#             contains no samples or all samples and thus are          #
#             discarded.                                               #
#                                                                      #
#    seed     An integer object.                                       #
#             If seed is provided as input.                            #
#                                                                      #
#----------------------------------------------------------------------#
subgroup_detect <- function(outcome,
                            propen,
                            data,
                            K = 1000L,
                            M = 1000L,
                            seed = NULL) {

  #------------------------------------------------------------------#
  # data must be a data.frame                                        #
  #------------------------------------------------------------------#
  if( !is(data, "data.frame") ) stop("data must be a data.frame.")

  #------------------------------------------------------------------#
  # data must be complete                                            #
  #------------------------------------------------------------------#
  if( any(is.na(data)) ) stop("data cannot contain missing values.")

  #------------------------------------------------------------------#
  # outcome must be a formula object.                                #
  #------------------------------------------------------------------#
  if( !is(outcome,"formula") ) stop("outcome must be a formula object.")
  toutcome <- terms(outcome)

  #------------------------------------------------------------------#
  # All variables of outcome must be provided in data.               #
  #------------------------------------------------------------------#
  mfo <- try(stats::model.frame(formula = toutcome, data = data))
  if( is(mfo, "try-error") ) {
    stop("Not all variables of outcome model are found in data.")
  }

  #------------------------------------------------------------------#
  # outcome must include an intercept.                               #
  #------------------------------------------------------------------#
  if( attr(toutcome,"intercept") < 0.5 ) {
    stop("outcome model must include an intercept.")
  }

  #------------------------------------------------------------------#
  # Response must be explicitly given on lhs of formula              #
  #------------------------------------------------------------------#
  y <- try(stats::model.response(data = mfo))
  if( is(y, "try-error") || is(y,"NULL") ) {
    stop("Response must be given on lhs of outcome model.")
  }

  #------------------------------------------------------------------#
  # Response must be continuous.                                     #
  #------------------------------------------------------------------#
  if( length(unique(y)) / length(y) < 0.10 ) {
    warning("Response variable may not be continuous.")
  }

  #------------------------------------------------------------------#
  # propen must be a formula object.                                 #
  #------------------------------------------------------------------#
  if( !is(propen,"formula") ) stop("propen must be a formula object.")
  tpropen <- terms(propen)

  #------------------------------------------------------------------#
  # All variables of propen model must be provided in data.          #
  #------------------------------------------------------------------#
  mfp <- try(stats::model.frame(formula = tpropen, data = data))
  if( is(mfp, "try-error") ) {
    stop("Not all variables of propen are found in data.")
  }

  #------------------------------------------------------------------#
  # Treatment must be explicitly given on lhs of formula             #
  #------------------------------------------------------------------#
  tx <- try(stats::model.response(data = mfp))
  if( is(tx, "try-error")  || is(tx,"NULL")) {
    stop("Treatment must be given on lhs of propen model.")
  }

  #------------------------------------------------------------------#
  # Treatment must be binary.                                        #
  #------------------------------------------------------------------#
  if( length(unique(round(tx,4L))) != 2L ) {
    stop("Treatment variable must be binary.")
  }

  #------------------------------------------------------------------#
  # M must be >= 1000                                                #
  #------------------------------------------------------------------#
  if( M < 999.5 ) {
    M <- 1000L
    warning("M reset to minimum value of 1000.")
  }

  #------------------------------------------------------------------#
  # If a seed is specified, set seed.                                #
  #------------------------------------------------------------------#
  if( !is(seed,"NULL") )  {
    if( !is(seed,"integer") ) seed <- as.integer(seed)
    set.seed(seed)
  }

  #------------------------------------------------------------------#
  # parameter estimation for outcome regression                      #
  #------------------------------------------------------------------#
  model_mu <- try(stats::lm(formula = outcome, data = data))
  if( is(model_mu, "try-error") ) {
    stop("Unable to obtain parameter estimates for outcome using lm.")
  }
  beta_est <- stats::coef(model_mu)
  if( any(is.na(beta_est)) ) {
    stop("Some outcome model parameters returned as NA.")
  }
  mu_fit <- stats::fitted.values(model_mu)
  XO <- stats::model.matrix(model_mu)

  #------------------------------------------------------------------#
  # Warn if K is smaller than recommended value.                     #
  #------------------------------------------------------------------#
  p <- ncol(XO)
  n <- nrow(XO)

  if( K < 10^p ) {
    warning("K is smaller than recommended. Minimum ~ 10^p.")
  }

  #------------------------------------------------------------------#
  # parameter estimation for propensity                              #
  #------------------------------------------------------------------#
  model_pi <- try(stats::glm(propen, 
                             family = binomial(link="logit"),
                             data = data))
  if( is(model_pi, "try-error") ) {
    stop("Unable to obtain parameter estimates for propensity using glm.")
  }
  gamma_est <- stats::coef(model_pi)
  if( any(is.na(gamma_est)) ) {
    stop("Some outcome model parameters returned as NA.")
  }
  pi_fit <- stats::fitted.values(model_pi)
  XP <- stats::model.matrix(model_pi)

  #------------------------------------------------------------------#
  # influence functions estimate                                     #
  #------------------------------------------------------------------#
  spsi2 <- t(XO * {y - mu_fit})
  spsi3 <- t(XP * {tx - pi_fit})

  #------------------------------------------------------------------#
  # derivatives of influence functions                               #
  #------------------------------------------------------------------#
  dpsi22 <- -crossprod(XO, XO) / n
  dpi <- pi_fit^2 * exp(-drop(XP %*% gamma_est))
  dpsi33 <- -crossprod(XP * dpi, XP) / n

  #------------------------------------------------------------------#
  # influence function of test statistic                             #
  #------------------------------------------------------------------#
  res1 <- try(solve(a = dpsi22, b = spsi2))
  if( is(res1, "try-error") ) {
    stop("Unable to invert matrix derivative of influence function.")
  }

  res2 <- try(solve(a = dpsi33, b = spsi3))
  if( is(res2, "try-error") ) {
    stop("Unable to invert matrix derivative of influence function.")
  }

  #------------------------------------------------------------------#
  # Initialize storage variables                                     #
  #------------------------------------------------------------------#
  teststat <- -Inf
  sup_ts <- 0L

  teststat_p <- matrix(data = NA, nrow = M, ncol = K)

  theta_all <- matrix(data = NA, 
                      nrow = K, 
                      ncol = p, 
                      dimnames = list(NULL, colnames(XO)))

  prop <- 0L

  G <- matrix(data = rnorm(M*n), nrow = n, ncol = M)

  for( i in 1L:K ) {

    #--------------------------------------------------------------#
    # uniformly distributed on unit ball                           #
    #--------------------------------------------------------------#
    temp <- stats::rnorm(n = p, mean = 0.0, sd = 1.0)
    theta <- temp / sqrt(sum(temp^2)) 
    theta_all[i,] <- theta

    #--------------------------------------------------------------#
    # Calculate indicator function                                 #
    #--------------------------------------------------------------#
    subgroup <- drop(XO %*% theta) >= 0.0
    if( !any(subgroup) | all(subgroup) ) next

    prop <- prop + 1L

    #--------------------------------------------------------------#
    # calculate test statistic                                     #
    #--------------------------------------------------------------#
    ta <- tx - pi_fit
    ty <- y - mu_fit

    spsi1 <- ta * ty * subgroup
    dpsi12 <- -colSums(XO * subgroup * ta) / n
    dpsi13 <- -colSums(XP * subgroup * ty * dpi) / n

    psi_t <- spsi1 - drop(dpsi12 %*% res1) - drop(dpsi13 %*% res2)

    #--------------------------------------------------------------#
    # variance estimation of test statistic                        #
    #--------------------------------------------------------------#
    V_t <- 1.0 / sum(psi_t^2)

    #--------------------------------------------------------------#
    # standardized test statistic                                  #
    #--------------------------------------------------------------#
    temp <- {sum(spsi1)}^2 * V_t
    if( temp > teststat ) {
      teststat <- temp
      sup_ts <- i
    }

    #--------------------------------------------------------------#
    # empirical distribution of the test statistic                 #
    #--------------------------------------------------------------#
    teststat_p[,i] <- {colSums(psi_t * G)}^2 * V_t

  }

  #------------------------------------------------------------------#
  # distribution of test statistic                                   #
  #------------------------------------------------------------------#
  dissup_ts <- apply(X = teststat_p,
                     MARGIN = 1L,
                     FUN = max,
                     na.rm = TRUE)

  #------------------------------------------------------------------#
  # p-value                                                          #
  #------------------------------------------------------------------#
  p_value <- sum(dissup_ts >= teststat) / length(dissup_ts)

  #------------------------------------------------------------------#
  # change plane parameter estimate                                  #
  #------------------------------------------------------------------#
  theta <- theta_all[sup_ts,]

  #------------------------------------------------------------------#
  # output                                                           #
  #------------------------------------------------------------------#
  obj <- list("outcome" = model_mu,
              "propen" = model_pi,
              "p_value" = p_value,
              "theta" = theta,
              "prop" = prop)

  if( !is(seed, "NULL") ) obj$seed <- seed

  return(obj)
}

