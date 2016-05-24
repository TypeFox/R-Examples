
trunc01 <- function(x){

  pmax(pmin(x, 1), 0)

}

expit <- function(x) exp(x)/(1 + exp(x))

#' Generate sample data used for testing
#'
#' @param n Integer, the sample size
#'
#' @export

generate_example_data <- function(n){

  Z <- rbinom(n, 1, .5)
  X <- rnorm(n)  ## pre-treatment covariate

  ## generate S conditional on Z

  S.0 <- X + rnorm(n, sd = .1)
  S.1 <- 1 + X + rnorm(n, sd = .1)

  beta <- .75

  risk.obs <- (-1 - 0.0 * S.1 - 0 * Z - beta * S.1 * Z)
  risk.0 <- (-1 - 0.0 * S.1)
  risk.1 <- (-1 - beta * S.1)

  time.0 <- rexp(n, 1/exp(risk.0))
  time.1 <- rexp(n, 1/exp(risk.1))
  time.obs <- ifelse(Z == 1, time.1, time.0)

  event.0 <- rbinom(n, 1, .8)
  event.1 <- rbinom(n, 1, .8)

  event.0[time.0 > 150] <- 0
  event.1[time.1 > 150] <- 0
  event.obs <- ifelse(Z == 1, event.1, event.0)

  Y.0 <- rbinom(n, 1, expit(risk.0))
  Y.1 <- rbinom(n, 1, expit(risk.1))
  Y.obs <- ifelse(Z == 1, Y.1, Y.0)

  ## CPV measure noisy S.1 at the end of the study for placebo subjects non-event

  CPV <- S.1 + rnorm(n, sd = .1)
  CPV[Z == 1 | Y.obs == 1] <- NA

  ## BSM measure noisy S.0 at start

  BSM <- S.0 + rnorm(n, sd = .1)
  S.1[Z == 0] <- NA
  S.0[Z == 1] <- NA
  S.obs <- ifelse(Z == 1, S.1, S.0)

  ## make categorical variables
  qwantz <- c(-Inf, quantile(c(S.0, S.1), c(.25, .5, .75), na.rm = TRUE), Inf)
  S.1.cat <- cut(S.1, qwantz)
  S.0.cat <- cut(S.0, qwantz)
  S.obs.cat <- cut(S.obs, qwantz)

  BIP.cat <- cut(X, c(-Inf, quantile(X, c(.25, .5, .75)), Inf))

  data.frame(Z, BIP = X, CPV, BSM, S.obs, time.obs, event.obs, Y.obs, S.obs.cat,  BIP.cat)

}






