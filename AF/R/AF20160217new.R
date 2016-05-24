############## Common functions ##############
globalVariables(".SD")
is.binary <- function(v) {
  x <- unique(v)
  if((length(x) - sum(is.na(x)) == 2L) & (max(x) == 1 & min(x) == 0)) TRUE
  else FALSE
}

############## AF function for cross-sectional sampling design #####################
#' @title Attributable fraction for cross-sectional sampling designs.
#' @description \code{AF.cs} estimates the model-based adjusted attributable fraction for data from cross-sectional sampling designs.
#' @param formula an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model used for adjusting for confounders. The exposure and confounders should be specified as independent (right-hand side) variables. The outcome should be specified as dependent (left-hand side) variable. The formula is used to fit a logistic regression by \code{\link[stats]{glm}}.
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered.
#' @return \item{AF.est}{estimated attributable fraction.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta method with the sandwich formula.}
#' @return \item{P.est}{estimated factual proportion of cases; \eqn{Pr(Y=1)}.}
#' @return \item{P.var}{estimated variance of \code{P.est}. The variance is obtained by the sandwich formula.}
#' @return \item{P0.est}{estimated counterfactual proportion of cases if exposure would be eliminated; \eqn{Pr(Y_0=1)}{Pr(Y0=1)}.}
#' @return \item{P0.var}{estimated variance of \code{P0.est}. The variance is obtained by the sandwich formula.}
#' @return \item{fit}{the fitted model. Fitted using logistic regression, \code{\link{glm}}.}
#' @details \code{Af.cs} estimates the attributable fraction for a binary outcome \code{Y}
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by logistic regression (\code{\link{glm}}).
#' Let the AF be defined as
#' \deqn{AF=1-\frac{Pr(Y_0=1)}{Pr(Y=1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denotes the counterfactual probability of the outcome if
#' the exposure would have been eliminated from the population and \eqn{Pr(Y = 1)} denotes the factual probability of the outcome.
#' If \code{Z} is sufficient for confounding control, then \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be expressed as
#'  \eqn{E_Z\{Pr(Y=1\mid{X=0,Z})\}.}{E_z{Pr(Y = 1  |X = 0,Z)}.}
#' The function uses logistic regression to estimate \eqn{Pr(Y=1\mid{X=0,Z})}{Pr(Y=1|X=0,Z)}, and the marginal sample distribution of \code{Z}
#' to approximate the outer expectation (\enc{Sjölander}{Sjolander} and Vansteelandt, 2012).
#' If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link{glm}} used for fitting the logistic regression model.
#' @references Greenland, S. and Drescher, K. (1993). Maximum Likelihood Estimation of the Attributable Fraction from logistic Models. \emph{Biometrics} \bold{49}, 865-872.
#' @references \enc{Sjölander}{Sjolander}, A. and Vansteelandt, S. (2011). Doubly robust estimation of attributable fractions. \emph{Biostatistics} \bold{12}, 112-121.
#' @examples
#' # Simulate a cross-sectional sample
#' expit <- function(x) 1 / (1 + exp( - x))
#' n <- 1000
#' Z <- rnorm(n = n)
#' X <- rbinom(n = n, size = 1, prob = expit(Z))
#' Y <- rbinom(n = n, size = 1, prob = expit(Z + X))
#'
#' # Example 1: non clustered data from a cross-sectional sampling design
#' data <- data.frame(Y, X, Z)
#' AF.est.cs <- AF.cs(formula = Y ~ X + Z + X * Z, data = data, exposure = "X")
#' summary(AF.est.cs)
#'
#' # Example 2: clustered data from a cross-sectional sampling design
#' # Duplicate observations in order to create clustered data
#' id <- rep(1:n, 2)
#' data <- data.frame(id = id, Y = c(Y, Y), X = c(X, X), Z = c(Z, Z))
#' AF.est.cs.clust <- AF.cs(formula = Y ~ X + Z + X * Z, data = data,
#'                          exposure = "X", clusterid = "id")
#' summary(AF.est.cs.clust)
#' @importFrom stats aggregate ave binomial coef delete.response family glm model.matrix pnorm predict qnorm residuals stepfun terms var vcov
#' @export
AF.cs<- function(formula, data, exposure, clusterid){
  call <- match.call()
  mm <- match(c("formula", "data", "exposure", "clusterid"), names(call), 0L)
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  outcome <- as.character(terms(formula)[[2]])
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  if(missing(clusterid)) n.cluster <- 0
  else {
    n.cluster <- length(unique(data[, clusterid]))
    data <- data[order(data[, clusterid]), ]
  }
  ## Checks ##
  if(!is.binary(data[, outcome]))
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  ## Counterfactual dataset ##
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  fit <- glm(formula = formula, family = binomial, data = data)
  npar <- length(fit$coef)
  ## Design matrices ##
  design <- model.matrix(object = delete.response(terms(fit)), data = data)
  design0 <- model.matrix(object = delete.response(terms(fit)), data = data0)
  #### Meat: score equations ####
  ## Score equation 1 ##
  score.P <- data[, outcome]
  pred.Y  <- predict(fit, newdata = data, type = "response")
  ## Score equation 2 ##
  score.P0 <- predict(fit, newdata = data0, type = "response")
  ## Score equation 3 ##
  score.beta <- design * (score.P - pred.Y)
  ### Meat ###
  score.equations <- cbind(score.P, score.P0, score.beta)
  if (!missing(clusterid)){
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
  }
  meat <- var(score.equations, na.rm = TRUE)
  #### Bread: hessian of score equations ####
  ## Hessian of score equation 1 ##
  hessian.P <- matrix(c(- 1, 0, rep(0,npar)), nrow = 1, ncol = 2 + npar)
  ## Hessian of score equation 2 ##
  g <- family(fit)$mu.eta
  dmu.deta <- g(predict(object = fit, newdata = data0))
  deta.dbeta <- design0
  dmu.dbeta <- dmu.deta * deta.dbeta
  hessian.P0 <- matrix(c(0, - 1, colMeans(dmu.dbeta)), nrow = 1, ncol = 2 + npar)
  ## Hessian of score equation 3 ##
  hessian.beta <- cbind(matrix(rep(0, npar * 2), nrow = npar, ncol = 2)
                        , - solve(vcov(object = fit)) / n)
  ### Bread ###
  bread <- rbind(hessian.P, hessian.P0, hessian.beta)
  #### Sandwich ####
  if (!missing(clusterid))
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n^2 ) [1:2, 1:2]
  else
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
  #### Point estimate of AF ####
  P.est  <- mean(score.P, na.rm = TRUE)
  P0.est <- mean(score.P0, na.rm = TRUE)
  AF.est <- 1 - P0.est / P.est
  ## Delta method for variance estimate ##
  gradient <- as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
  AF.var <- t(gradient) %*% sandwich %*% gradient
  P.var <- sandwich[1, 1]
  P0.var <- sandwich[2, 2]
  #### Output ####
  out <- c(list(AF.est = AF.est, AF.var = AF.var, P.est = P.est, P0.est = P0.est, P.var = P.var,
                P0.var = P0.var, fitcall = fit$call, call = call, exposure = exposure, outcome = outcome,
                fit = fit, sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster))
  class(out) <- "AF"
  return(out)
}

############## AF function for cohort time-to-event outcomes #####################
#' @title Attributable fraction function for cohort sampling designs with time-to-event outcomes.
#' @description \code{AF.ch} estimates the model-based adjusted attributable fraction function for data from cohort sampling designs with time-to-event outcomes.
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function (\code{\link[survival]{Surv}}). The exposure and confounders should be specified as independent (right-hand side) variables. The time-to-event outcome should be specified by the survival object. The formula is used to fit a Cox proportional hazards model.
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param ties a character string specifying the method for tie handling. If there are no tied death times all the methods are equivalent. Uses the Breslow method by default.
#' @param times a scalar or vector of time points specified by the user for which the attributable fraction function is estimated. If not specified the observed death times will be used.
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered.
#' @return \item{AF.est}{estimated attributable fraction function for every time point specified by \code{times}.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta methods with the sandwich formula.}
#' @return \item{S.est}{estimated factual survival function; \eqn{S(t)}.}
#' @return \item{S.var}{estimated variance of \code{S.est}. The variance is obtained by the sandwich formula.}
#' @return \item{S0.est}{estimated counterfactual survival function if exposure would be eliminated; \eqn{S_0(t)}{S0(t)}.}
#' @return \item{S0.var}{estimated variance of \code{S0.est}. The variance is obtained by the sandwich formula.}
#' @return \item{fit}{the fitted model. Fitted using Cox proportional hazard, \code{\link[survival]{coxph}}.}
#' @details \code{Af.ch} estimates the attributable fraction for a time-to-event outcome
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population. The estimate is adjusted for confounders \code{Z}
#' by the Cox proportional hazards model (\code{\link[survival]{coxph}}). Let the AF function be defined as
#' \deqn{AF=1-\frac{\{1-S_0(t)\}}{\{1-S(t)\}}}{AF = 1 - {1 - S0(t)} / {1 - S(t)}}
#' where \eqn{S_0(t)}{S0(t)} denotes the counterfactual survival function for the event if
#' the exposure would have been eliminated from the population at baseline and \eqn{S(t)} denotes the factual survival function.
#' If \code{Z} is sufficient for confounding control, then \eqn{S_0(t)}{S0(t)} can be expressed as \eqn{E_Z\{S(t\mid{X=0,Z })\}}{E_z{S(t|X=0,Z)}}.
#' The function uses Cox proportional hazards regression to estimate \eqn{S(t\mid{X=0,Z})}{S(t|X=0,Z)}, and the marginal sample distribution of \code{Z}
#' to approximate the outer expectation (\enc{Sjölander}{Sjolander} and Vansteelandt, 2014).  If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link[survival]{coxph}} and \code{\link[survival]{Surv}} used for fitting the Cox proportional hazards model.
#' @references Chen, L., Lin, D. Y., and Zeng, D. (2010). Attributable fraction functions for censored event times. \emph{Biometrika} \bold{97}, 713-726.
#' @references \enc{Sjölander}{Sjolander}, A. and Vansteelandt, S. (2014). Doubly robust estimation of attributable fractions in survival analysis. \emph{Statistical Methods in Medical Research}. doi: 10.1177/0962280214564003.
#' @examples
#' # Simulate a sample from a cohort sampling design with time-to-event outcome
#' expit <- function(x) 1 / (1 + exp( - x))
#' n <- 500
#' time <- c(seq(from = 0.2, to = 1, by = 0.2))
#' Z <- rnorm(n = n)
#' X <- rbinom(n = n, size = 1, prob = expit(Z))
#' Tim <- rexp(n = n, rate = exp(X + Z))
#' C <- rexp(n = n, rate = exp(X + Z))
#' Tobs <- pmin(Tim, C)
#' D <- as.numeric(Tobs < C)
#' #Ties created by rounding
#' Tobs <- round(Tobs, digits = 2)
#'
#' # Example 1: non clustered data from a cohort sampling design with time-to-event outcomes
#' data <- data.frame(Tobs, D, X,  Z)
#' AF.est.ch <- AF.ch(formula = Surv(Tobs, D) ~ X + Z + X * Z, data = data,
#'                    exposure = "X", times = time)
#' summary(AF.est.ch)
#'
#' # Example 2: clustered data from a cohort sampling design with time-to-event outcomes
#' # Duplicate observations in order to create clustered data
#' id <- rep(1:n, 2)
#' data <- data.frame(Tobs = c(Tobs, Tobs), D = c(D, D), X = c(X, X), Z = c(Z, Z), id = id)
#' AF.est.ch.clust <- AF.ch(formula = Surv(Tobs, D) ~ X + Z + X * Z, data = data,
#'                          exposure = "X", times = time, clusterid = "id")
#' summary(AF.est.ch.clust)
#' plot(AF.est.ch.clust, CI = TRUE)
#' @import survival data.table
#' @export
AF.ch <- function(formula, data, exposure, ties="breslow",
                  times, clusterid){
  call <- match.call()
  mm <- match(c("formula", "data", "exposure", "ties", "times", "clusterid"), names(call), 0L)
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  ## If times is missing ##
  if(missing(times))
    times <- fit.detail$time
  ## Checks ##
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  if(missing(clusterid)) n.cluster <- 0
  else n.cluster <- length(unique(data[, clusterid]))
  ## Find names of end variable and event variable
  rr <- rownames(attr(terms(formula), "factors"))[1]
  temp <- gregexpr(", ", rr)[[1]]
  if(length(temp == 1)){
    endvar <- substr(rr, 6, temp[1] - 1)
    eventvar <- substr(rr, temp[1] + 2, nchar(rr) - 1)
  }
  if(length(temp) == 2){
    endvar <- substr(rr, temp[1] + 2, temp[2] - 1)
    eventvar <- substr(rr, temp[2] + 2, nchar(rr) - 1)
  }
  n <- nrow(data)
  n.cases <- sum(data[, eventvar])
  # Sort on "end-variable"
  data <- data[order(data[, endvar]), ]
  # Create dataset data0 for counterfactual X=0
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  ## Fit a Cox PH model ##
  environment(formula) <- new.env()
  fit <- coxph(formula = formula, data = data, ties = "breslow")
  npar <- length(fit$coef)
  fit.detail <- coxph.detail(object = fit)
  ## Design matrices ##
  design <- as.matrix(model.matrix(object = delete.response(terms(fit)), data = data)[, -1])
  design0 <- as.matrix(model.matrix(object = delete.response(terms(fit)), data = data0)[, -1])
  ### Estimate the survival functions ###
  ## Hazard increment ##
  dH0 <- fit.detail$hazard
  H0 <- cumsum(dH0)
  ## Baseline hazard function ##
  H0step <- stepfun(fit.detail$time, c(0, H0))
  H0res <- rep(0, n)
  dH0.untied <- rep(dH0, fit.detail$nevent) / rep(fit.detail$nevent, fit.detail$nevent)
  H0res[data[, eventvar] == 1] <- dH0.untied * n #handle ties
  #H0res[data[, eventvar] == 1] <- dH0 * n
  ## Predict based on the Cox PH model ##
  epred <- predict(object = fit, newdata = data, type = "risk")
  epred0 <- predict(object = fit, newdata = data0, type = "risk")
  ### Meat ###
  ## Score equation 4 ## for the Cox PH model (made outside of loop)
  score.beta <- residuals(object = fit, type = "score")
  ## Weighted mean of the variable at event for all at risk at that time ##
  E <- matrix(0, nrow = n, ncol = npar)
  means <- as.matrix(fit.detail$means)
  means <- means[rep(1:nrow(means), fit.detail$nevent), ] #handle ties
  E[data[, eventvar] == 1, ] <- means
  #E[data[, eventvar] == 1, ] <- fit.detail$means
  ## One point and variance estimate for each time t in times ##
  S.est <- vector(length = length(times))
  S0.est <- vector(length = length(times))
  AF.var <- vector(length = length(times))
  S.var <- vector(length = length(times))
  S0.var <- vector(length = length(times))
  
  # Loop over all t in times
  for (i in 1:length(times)){
    t <- times[i]
    #### Meat: score equations ####
    ## Score equation 1 ## for the factual survival function
    score.S <- exp( - H0step(t) * epred)
    ## Score equation 2 ## for the counterfactual survival function
    score.S0 <- exp( - H0step(t) * epred0)
    ## Score equation 3 ##  for the Breslow estimator
    score.H0 <- H0res * (data[, endvar] <= t)
    ## Score equation 4 ## for the Cox PH model (made outside of loop)
    ### Meat ###
    score.equations <- cbind(score.S, score.S0, score.H0, score.beta)
    if (!missing(clusterid)){
      #score.equations <-aggregate(score.equations, by = list(data[, clusterid]), sum)[, - 1]
      score.equations <- data.table(score.equations)
      score.equations <- score.equations[, j=lapply(.SD,sum), by=clusterid]
      score.equations <- as.matrix(score.equations)
      score.equations <- score.equations[, -1]
    }
    
    meat <- var(score.equations, na.rm = TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    hessian.S <- c(-1, 0, mean(epred * score.S), colMeans(design * H0step(t) * epred * score.S))
    ## Hessian of score equation 2 ##
    hessian.S0 <- c(0, -1, mean(epred0 * score.S0), colMeans(design0 * H0step(t) * epred0 * score.S0))
    ## Hessian of score equation 3 ##
    hessian.H0 <- c(rep(0,2), - 1, - colMeans(E * score.H0, na.rm = TRUE))
    ## Hessian of score equation 4 ##
    hessian.beta <- cbind(matrix(0, nrow = npar, ncol = 3), - solve(vcov(object = fit)) / n)
    ### Bread ###
    bread<-rbind(hessian.S, hessian.S0, hessian.H0, hessian.beta)
    ### Sandwich ###
    if (!missing(clusterid))
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n^2 ) [1:2, 1:2]
    else
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    #### For point estimate ####
    S.est[i] <- mean(x = score.S, na.rm = TRUE)
    S0.est[i] <- mean(x = score.S0, na.rm = TRUE)
    #### Estimate of variance using the delta method ####
    gradient <- as.matrix(c( - (1 - S0.est[i]) / (1 - S.est[i]) ^ 2, 1 / (1 - S.est[i]))
                          , nrow = 2, ncol = 1)
    AF.var[i] <- t(gradient) %*% sandwich %*% gradient
    S.var[i] <- sandwich[1, 1]
    S0.var[i] <- sandwich[2, 2]
  }
  ### The AF function estimate ###
  AF.est <- 1 - (1 - S0.est) / (1 - S.est)
  #### Output ####
  #func <- AF.cc
  out <- c(list(AF.est = AF.est, AF.var = AF.var, S.est = S.est,
                S0.est = S0.est, S.var = S.var, S0.var = S0.var,
                fitcall = fit$call, call = call, exposure = exposure, outcome = eventvar, fit = fit,
                sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster,  times = times))
  class(out) <- "AF"
  return(out)
}

############## AF function for matched and unmatched case-control #####################
#' @title Attributable fraction for mached and non-matched case-control sampling designs.
#' @description \code{AF.cc} estimates the model-based adjusted attributable fraction for data from matched and non-matched case-control sampling designs.
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model used for confounder adjustment. The exposure and confounders should be specified as independent (right-hand side) variables. The outcome should be specified as dependent (left-hand side) variable. The formula is used to fit a logistic regression by \code{\link[stats]{glm}} for non-matched case-control and conditional logistic regression by \code{\link[drgee]{gee}} (in package \code{\link[drgee]{drgee}}) for matched case-control.
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param matched a logical that specifies if the sampling design is matched (TRUE) or non-matched (FALSE) case-control. Default setting is non-matched (\code{matched = FALSE}).
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered (e.g. matched).
#' @return \item{AF.est}{estimated attributable fraction.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta methods with the sandwich formula.}
#' @return \item{log.or}{a vector of the estimated log odds ratio for every individual. \code{log.or} contains the estimated coefficient for the exposure variable \code{X} for every level of the confounder \code{Z} as specified by the user in the formula. If the model to be estimated is
#'  \deqn{logit\{Pr(Y=1|X,Z)\} = \alpha+\beta{X}+\gamma{Z}}{logit {Pr(Y=1|X,Z)} = \alpha + \beta X + \gamma Z}
#'   then \code{log.or} is the estimate of \eqn{\beta}.
#'   If the model to be estimated is
#'   \deqn{logit\{Pr(Y=1|X,Z)\}=\alpha+\beta{X}+\gamma{Z}+\psi{XZ}}{logit{Pr(Y=1|X,Z)} = \alpha + \beta X +\gamma Z +\psi XZ}
#'   then \code{log.odds} is the estimate of
#'    \eqn{\beta + \psi{Z}}{\beta + \psi Z}.}
#' @return \item{fit}{the fitted model. Fitted using logistic regression, \code{\link{glm}}, for non-matched case-control and conditional logistic regression, \code{\link[drgee]{gee}}, for matched case-control.}
#' @details \code{Af.cc} estimates the attributable fraction for a binary outcome \code{Y}
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by logistic regression for unmatched case-control (\code{\link[stats]{glm}}) and conditional logistic regression for matched case-control (\code{\link[drgee]{gee}}).
#' The estimation assumes that the outcome is rare so that the risk ratio can be approximated by the odds ratio, for details see Bruzzi et. al.
#' Let the AF be defined as
#' \deqn{AF = 1 - \frac{Pr(Y_0=1)}{Pr(Y = 1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denotes the counterfactual probability of the outcome if
#' the exposure would have been eliminated from the population. If \code{Z} is sufficient for confounding control then the probability \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be expressed as
#' \deqn{Pr(Y_0=1)=E_Z\{Pr(Y=1\mid{X}=0,Z)\}.}{Pr(Y0=1) = E_z{Pr(Y = 1 | X = 0, Z)}.}
#' Using Bayes' theorem this implies that the AF can be expressed as
#' \deqn{AF = 1-\frac{E_Z\{Pr(Y=1\mid X=0,Z)\}}{Pr(Y=1)}=1-E_Z\{RR^{-X}(Z)\mid{Y = 1}\}}{
#' AF = 1 - E_z{Pr( Y = 1 | X = 0, Z)} / Pr(Y = 1) = 1 - E_z{RR^{-X} (Z) | Y = 1}}
#' where \eqn{RR(Z)} is the risk ratio \deqn{\frac{Pr(Y=1\mid{X=1,Z})}{Pr(Y=1\mid{X=0,Z})}.}{Pr(Y = 1 | X = 1,Z)/Pr(Y=1 | X = 0, Z).}
#' Moreover, the risk ratio can be approximated by the odds ratio if the outcome is rare. Thus,
#' \deqn{ AF \approx 1 - E_Z\{OR^{-X}(Z)\mid{Y = 1}\}.}{AF is approximately 1 - E_z{OR^{-X}(Z) | Y = 1}.}
#' The odds ratio is estimated by logistic regression or conditional logistic regression.
#' If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link[stats]{glm}} and \code{\link[drgee]{gee}} used for fitting the logistic regression model (for non-matched case-control) and the conditional logistic regression model (for matched case-control).
#' @references Bruzzi, P., Green, S. B., Byar, D., Brinton, L. A., and Schairer, C. (1985). Estimating the population attributable risk for multiple risk factors using case-control data. \emph{American Journal of Epidemiology} \bold{122}, 904-914.
#' @examples
#' expit <- function(x) 1 / (1 + exp( - x))
#' NN <- 1000000
#' n <- 500
#'
#' # Example 1: non matched case-control
#' # Simulate a sample from a non matched case-control sampling design
#' # Make the outcome a rare event by setting the intercept to -6
#' intercept <- -6
#' Z <- rnorm(n = NN)
#' X <- rbinom(n = NN, size = 1, prob = expit(Z))
#' Y <- rbinom(n = NN, size = 1, prob = expit(intercept + X + Z))
#' population <- data.frame(Z, X, Y)
#' Case <- which(population$Y == 1)
#' Control <- which(population$Y == 0)
#' # Sample cases and controls from the population
#' case <- sample(Case, n)
#' control <- sample(Control, n)
#' data <- population[c(case, control), ]
#' AF.est.cc <- AF.cc(formula = Y ~ X + Z + X * Z, data = data, exposure = "X")
#' summary(AF.est.cc)
#'
#' # Example 2: matched case-control
#' # Duplicate observations in order to create a matched data sample
#' # Create an unobserved confounder U common for each pair of individuals
#' U  <- rnorm(n = NN)
#' Z1 <- rnorm(n = NN)
#' Z2 <- rnorm(n = NN)
#' X1 <- rbinom(n = NN, size = 1, prob = expit(U + Z1))
#' X2 <- rbinom(n = NN, size = 1, prob = expit(U + Z2))
#' Y1 <- rbinom(n = NN, size = 1, prob = expit(intercept + U + Z1 + X1))
#' Y2 <- rbinom(n = NN, size = 1, prob = expit(intercept + U + Z2 + X2))
#' # Select discordant pairs
#' discordant <- which(Y1!=Y2)
#' id <- rep(1:n, 2)
#' # Sample from discordant pairs
#' incl <- sample(x = discordant, size = n, replace = TRUE)
#' data <- data.frame(id = id, Y = c(Y1[incl], Y2[incl]), X = c(X1[incl], X2[incl]),
#'                    Z = c(Z1[incl], Z2[incl]))
#' AF.est.cc.match <- AF.cc(formula = Y ~ X + Z + X * Z, data = data,
#'                          exposure = "X", clusterid = "id", matched = TRUE)
#' summary(AF.est.cc.match)
#' @import drgee
#' @export
AF.cc<-function(formula, data, exposure, clusterid,
                matched = FALSE){
  call <- match.call()
  mm <- match(c("formula", "data", "exposure", "clusterid", "matched"), names(call), 0L)
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  outcome <- as.character(terms(formula)[[2]])
  if(matched == TRUE){
    ni.vals <- ave(as.vector(data[, outcome]), data[, clusterid], FUN = function(y) {
      length(unique(y[which(!is.na(y))]))
    })
    compl.rows <- (ni.vals > 1)
    data <- data[compl.rows, ]
  }
  ## Checks ##
  if(is.binary(data[, outcome]) == FALSE)
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(is.binary(data[, exposure]) == FALSE)
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  if(missing(clusterid)) n.cluster <- 0
  else n.cluster <- length(unique(data[, clusterid]))
  #### Methods for non-matched or matched sampling designs ####
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  if (!missing(clusterid))
    data <- data[order(data[,clusterid]), ]
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  if(matched == FALSE)
    fit <- glm(formula = formula, family = binomial, data = data)
  if(matched == TRUE)
    fit <- gee(formula, link = "logit", data, cond = TRUE, clusterid = clusterid)
  npar <- length(fit$coef)
  ## Design matrices ##
  if(matched == FALSE){
    design <- as.matrix(model.matrix(object = delete.response(terms(fit)), data = data))
    design0 <- as.matrix(model.matrix(object = delete.response(terms(fit)), data = data0))
  }
  if(matched == TRUE){
    design <- as.matrix(model.matrix(object = formula, data = data)[, - 1])
    design0 <- as.matrix(model.matrix(object = formula, data = data0)[, - 1])
  }
  ## Create linear predictors to estimate the log odds ratio ##
  diff.design <- design0 - design
  linearpredictor <- design  %*% coef(fit)
  linearpredictor0 <- design0 %*% coef(fit)
  #log odds ratio#
  log.or <- linearpredictor - linearpredictor0
  ## Estimate approximate AF ##
  AF.est   <- 1 - sum(data[, outcome] * exp( - log.or)) / sum(data[, outcome])
  #### Meat: score equations ####
  ## Score equation 1 ## individual estimating equations of the estimate of AF
  score.AF <- data[, outcome] * (exp( - log.or) - AF.est)
  ## Score equation 2 ## individual estimating equations from conditional logistic reg.
  if(matched == FALSE)
    pred.diff <- data[, outcome] - predict(fit, newdata = data, type = "response")
  if(matched == TRUE)
    pred.diff <- fit$res
  score.beta <- design * pred.diff
  score.equations <- cbind(score.AF, score.beta)
  if (!missing(clusterid))
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
  meat <- var(score.equations, na.rm=TRUE)
  #### Bread: hessian of score equations ####
  ## Hessian of score equation 1 ##
  #### Estimating variance using Sandwich estimator ####
  hessian.AF1 <- - data[, outcome]
  
  
  hessian.AF2 <- (design0 - design) * as.vector(data[, outcome] * exp( - log.or))
  if (!missing(clusterid)){
    if(length(all.vars(formula[[3]]))>1){
      hessian.AF <- cbind(mean(aggregate(hessian.AF1, list(data[, clusterid]), sum)[, - 1], na.rm=TRUE)
                          , t(colMeans(aggregate(hessian.AF2
                                                 , list(data[, clusterid]), sum)[, - 1], na.rm = TRUE)))
    }
    if(length(all.vars(formula[[3]]))==1){
      hessian.AF <- cbind(mean(aggregate(hessian.AF1, list(data[, clusterid]), sum)[, - 1], na.rm=TRUE)
                          , t(mean(aggregate(hessian.AF2
                                             , list(data[, clusterid]), sum)[, - 1], na.rm = TRUE)))
    }
  }
  else
    hessian.AF <- cbind(mean(hessian.AF1), t(colMeans(hessian.AF2, na.rm = TRUE)))
  hessian.beta <- cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), - solve(vcov(object = fit)) / n)
  ### Bread ###
  bread <- rbind(hessian.AF, hessian.beta)
  #### Sandwich ####
  if (!missing(clusterid))
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n ^ 2 ) [1:2, 1:2]
  else
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
  AF.var <- sandwich[1, 1]
  #### Output ####
  out <- c(list(AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                fitcall = fit$call, call = call, exposure = exposure, outcome = outcome, fit = fit,
                sandwich = sandwich, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster))
  class(out) <- "AF"
  return(out)
}

############## Summary and print functions ##############
#' @export
print.AF<-function(x, ...){
  if(!x$n.cluster == 0) {
    Std.Error <- "Robust SE"
    se <- "cluster-robust standard error"
  }
  else {
    Std.Error <- "Std.Error"
    se <- "standard error"
  }
  cat("\nEstimated attributable fraction (AF) and", se, ":", "\n")
  cat("\n")
  table.est <- cbind(x$AF.est, sqrt(x$AF.var))
  colnames(table.est) <- c("AF", Std.Error)
  r <- rep("", , length(x$AF.est))
  rownames(table.est) <- c(r)
  modelcall <- as.character(x$fit$call[1])
  if(modelcall == "coxph") {
    table.est <- cbind(x$times, table.est)
    colnames(table.est) <- c("Time", "AF", Std.Error)
    print.default(table.est)
  }
  else {
    print.default(table.est)
  }
}

CI <- function(AF, Std.Error, confidence.level, CI.transform){
  if(CI.transform == "untransformed"){
    lower <- AF - abs(qnorm((1 - confidence.level) / 2)) * Std.Error
    upper <- AF + abs(qnorm((1 - confidence.level) / 2)) * Std.Error
  }
  if(CI.transform == "log"){
    lower <- AF * exp( - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / AF)
    upper <- AF * exp(abs(qnorm((1 - confidence.level) / 2)) * Std.Error / AF)
  }
  if(CI.transform == "logit"){
    logit <- function(x) log(x / (1 - x))
    lower <- exp(logit(AF) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))) / (1 + exp(logit(AF) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))))
    upper <- exp(logit(AF) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))) / (1 + exp(logit(AF) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))))
  }
  CI <- cbind(lower, upper)
  return(CI)
}

#' @title Summary function for objects of class "\code{AF}".
#' @description Gives a summary of the AF estimate(s) including z-value, p-value and confidence interval(s).
#' @param object an object of class \code{AF} from \code{\link{AF.cs}}, \code{\link{AF.ch}} or \code{\link{AF.cc}} functions.
#' @param confidence.level user-specified confidence level for the confidence intervals. If not specified it defaults to 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @param CI.transform user-specified transformation of the Wald confidence interval(s). Options are \code{untransformed}, \code{log} and \code{logit}. If not specified untransformed will be calculated.
#' @param digits maximum number of digits.
#' @param ... further arguments to be passed to the summary function. See \code{\link[base]{summary}}.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @export
summary.AF <- function(object, digits = max(3L, getOption("digits") - 3L),
                       confidence.level, CI.transform, ...){
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  se <- sqrt(object$AF.var)
  zvalue <- object$AF.est / sqrt(object$AF.var)
  pvalue <- 2 * pnorm( - abs(zvalue))
  confidence.interval <- CI(AF = object$AF.est, Std.Error = sqrt(object$AF.var),
                            confidence.level = confidence.level,
                            CI.transform = CI.transform)
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")
  
  if(!object$n.cluster == 0) Std.Error <- "Robust SE"
  else Std.Error <- "Std.Error"
  AF <- cbind(object$AF.est, se, zvalue, pvalue)
  colnames(AF) <- c("AF estimate", Std.Error, "z value", "Pr(>|z|)")
  
  modelcall <- as.character(object$fit$call[1])
  if(modelcall == "glm") method = "Logistic regression"
  if(modelcall == "coxph") method = "Cox Proportional Hazards model"
  if(modelcall == "gee") method = "Conditional logistic regression"
  
  fit <- summary(object$fit)
  
  if(modelcall == "coxph"){
    ans <- list(AF = AF, times = object$times,
                CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n,
                n.cases = object$n.cases, n.cluster = object$n.cluster,
                modelcall = modelcall, call = object$call, method = method, formula = object$formula,
                exposure = object$exposure, outcome = object$outcome, fit = fit,
                sandwich = object$sandwich)
  }
  else{
    ans <- list(AF = AF, CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n, n.cases = object$n.cases,
                n.cluster = object$n.cluster, modelcall = modelcall, call = object$call, method = method,
                formula = object$formula, exposure = object$exposure, outcome = object$outcome,
                fit = fit, sandwich = object$sandwich)
  }
  class(ans) <- "summary.AF"
  return(ans)
}

#' @export
print.summary.AF <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...){
  cat("Call: ", "\n")
  print.default(x$call)
  if(!x$n.cluster == 0) Std.Error <- "Robust SE"
  else Std.Error <- "Std.Error"
  if(x$CI.transform == "log") x$CI.transform <- "log transformed"
  if(x$CI.transform == "logit") x$CI.transform <- "logit transformed"
  level <- x$confidence.level * 100
  CI.text <- paste0(as.character(level),"%")
  cat("\nEstimated attributable fraction (AF) and", x$CI.transform, CI.text,  "Wald CI:", "\n")
  cat("\n")
  table.est <- cbind(x$AF, x$confidence.interval)
  colnames(table.est) <- c("AF", Std.Error, "z value", "Pr(>|z|)",
                           "Lower limit", "Upper limit")
  r <- rep("", , nrow(x$AF))
  rownames(table.est) <- c(r)
  modelcall <- as.character(x$fit$call[1])
  if(x$modelcall == "coxph") {
    table.est <- cbind(x$times, table.est)
    colnames(table.est) <- c("Time", "AF", Std.Error, "z value", "Pr(>|z|)",
                             "Lower limit", "Upper limit")
    print.default(table.est)
  }
  else {
    print.default(table.est)
  }
  cat("\nExposure", ":", x$exposure, "\n")
  
  if(x$modelcall == "coxph") outcome <- "Event   "
  else outcome <- "Outcome "
  #cat("\n")
  cat(outcome, ":", x$outcome, "\n")
  
  cat("\n")
  table.nr <- cbind(x$n.obs, x$n.cases)
  rownames(table.nr) <- c("")
  if(x$modelcall == "coxph") number <- "Events"
  else number <- "Cases"
  colnames(table.nr) <- c("Observations", number)
  if (x$n.cluster == 0) print.default(table.nr)
  else{
    table.nr.cluster <- cbind(table.nr, x$n.cluster)
    colnames(table.nr.cluster) <- c("Observations", number, "Clusters")
    print.default(table.nr.cluster)
  }
  cat("\nMethod for confounder adjustment: ", x$method, "\n")
  formula <- as.character(x$formula)
  cat("\nFormula: ", formula[2], formula[1], formula[3], "\n")
}

#' @title Plot function for objects of class "\code{AF}" from the function \code{\link{AF.ch}}.
#' @description Creates a simple scatterplot for the AF function with time sequence (specified by the user as \code{times} in the \code{\link{AF.ch}} function) on the x-axis and the AF function estimate on the y-axis.
#' @param x an object of class \code{AF} from the \code{\link{AF.ch}} function.
#' @param CI if TRUE confidence intervals are estimated and ploted in the graph.
#' @param confidence.level user-specified confidence level for the confidence intervals. If not specified it defaults to 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @param CI.transform user-specified transformation of the Wald confidence interval(s). Options are \code{untransformed}, \code{log} and \code{logit}. If not specified untransformed will be calculated.
#' @param xlab label on the x-axis. If not specified the label \emph{"Time"} will be displayed.
#' @param main main title of the plot. If not specified the lable \emph{"Estimate of the attributable fraction function"} will be displayed.
#' @param ... further arguments to be passed to the plot function. See \code{\link[graphics]{plot}}.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @importFrom graphics legend lines plot.default
#' @export
plot.AF <- function(x, CI = FALSE, confidence.level,
                    CI.transform, xlab, main, ...){
  modelcall <- as.character(x$fit$call[1])
  if(!modelcall == "coxph")
    stop("Plot function is only available for the attributable fraction function. That is objects from the AF.ch function", call. = FALSE)
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  if(missing(xlab)) xlab <- "Time"
  if(missing(main)) main <- "Estimate of the attributable fraction function"
  plot.default(x$times, x$AF.est, main = main,
               ylab = "Attributable fraction function" , xlab = xlab, pch = 19,
               lty = 1, type = "o", ...)
  if(CI == TRUE){
    confidence.interval <- CI(AF = x$AF.est, Std.Error = sqrt(x$AF.var),
                              confidence.level = confidence.level,
                              CI.transform = CI.transform)
    lines( x$times, confidence.interval[, 2], lty = 2)
    lines( x$times, confidence.interval[, 1], lty = 2)
    level <- confidence.level * 100
    CI <- paste0(as.character(level),"% Conf. Interval")
    if(CI.transform == "log") transform <- "(log transformed)"
    if(CI.transform == "logit") transform <- "(logit transformed)"
    if(CI.transform == "untransformed") transform <- "(untransformed)"
    legend("topright", legend = c("AF estimate", CI, transform), pch = c(19, NA, NA), lty = c(1, 2, 0), bty = "n")
  }
}

