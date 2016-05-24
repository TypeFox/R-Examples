#' Integration models
#'
#' Add integration model to a psdesign object
#'
#' @details This is a list of the available integration models. The fundamental problem in surrogate evaluation is that there are unobserved values of the counterfactual surrogate reponses S(1). In the estimated maximum likelihood framework, for subjects missing the S(1) values, we use an auxiliary pre-treatment variable or set of variables W that is observed for every subject to estimate the distribution of S(1) | W. Typically, this W is a BIP. Then for each missing S(1), we integrate likelihood contributions over each non-missing S(1) given their value of W, and average over the contributions.
#'
#' \itemize{
#' \item \link{integrate_parametric} This is a parametric integration model that fits a linear model for the mean of S(1) | W and assumes a Gaussian distribution.
#' \item \link{integrate_bivnorm} This is another parametric integration model that assumes that S(1) and W are jointly normally distributed. The user must specify their mean, variances and correlation.
#' \item \link{integrate_nonparametric} This is a non-parametric integration model that is only valid for categorical S(1) and W. It uses the observed proportions to estimate the joint distribution of S(1), W.
#' \item \link{integrate_semiparametric} This is a semi-parametric model that uses the semi-parametric location scale model of Heagerty and Pepe (1999). Models are specified for the location of S(1) | W and the scale of S(1) | W. Then integrations are drawn from the empirical distribution of the residuals from that model, which are then transformed to the appropriate location and scale.
#' }
#'
#'
#' @param psdesign A psdesign object
#' @param integration An integration object
#'
#' @export
#'
#' @examples
#'
#' test <- psdesign(generate_example_data(n = 100), Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
#' add_integration(test, integrate_parametric(S.1 ~ BIP))
#' test + integrate_parametric(S.1 ~ BIP)  # same as above
#'

add_integration <- function(psdesign, integration){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(integration, "integration"))

  psdesign <- integration(psdesign)
  psdesign

}

#' Add risk model to a psdesign object
#'
#' @details The risk model component specifies the likelihood for the data. This involves specifying the distribution of the outcome variable, whether it is binary or time-to-event, and specifying how the surrogate S(1) and the treatment Z interact and affect the outcome. We use the formula notation to be consistent with other regression type models in R. Below is a list of available risk models.
#'
#' \itemize{
#' \item \link{risk_binary} This is a generic risk model for binary outcomes. The user can specify the formula, and link function using either \link{risk.logit} for the logistic link, or \link{risk.probit} for the probit link. Custom link functions may also be specified, which take a single numeric vector argument, and returns a vector of corresponding probabilities.
#' \item \link{risk_weibull} This is a parameterization of the Weibull model for time-to-event outcomes that is consistent with that of \link{rweibull}. The user specifies the formula for the linear predictor of the scale parameter.
#' \item \link{risk_exponential} This is a simple exponential model for a time-to-event outcome.
#' \item \link{risk_poisson} This is a Poisson model for count outcomes. It allows for offsets in the formula.
#' }
#'
#' @param psdesign A psdesign object
#' @param riskmodel A risk model object, from the list above
#'
#' @export
#'
#' @examples
#' test <- psdesign(generate_example_data(n = 100), Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
#'      integrate_parametric(S.1 ~ BIP)
#' add_riskmodel(test, risk_binary())
#' test + risk_binary() # same as above

add_riskmodel <- function(psdesign, riskmodel){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(riskmodel, "riskmodel"))

  psdesign <- riskmodel(psdesign)
  psdesign

}

#' Estimate parameters
#'
#' @param psdesign A psdesign object, it must have risk model and integration model components
#' @param estimate An estimate object created by \link{ps_estimate}
#'
#' @examples
#' test <- psdesign(generate_example_data(n = 100), Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
#' test + integrate_parametric(S.1 ~ BIP) + risk_binary(D = 50) + ps_estimate(method = "BFGS")

add_estimate <- function(psdesign, estimate){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(estimate, "estimate"))

  psdesign <- estimate(psdesign)
  psdesign

}

#' Bootstrap resampling parameters
#'
#' @param psdesign A psdesign object, it must have risk model and integration model components
#' @param bootstrap A bootstrap object created by \link{ps_bootstrap}
#'
#' @examples
#' \dontrun{
#' test <- psdesign(generate_example_data(n = 100), Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
#' est1 <- test + integrate_parametric(S.1 ~ BIP) + risk_binary() + ps_estimate(method = "BFGS")
#' est1 + ps_bootstrap(method = "BFGS", start = est1$estimates$par, n.boots = 50)
#' }

add_bootstrap <- function(psdesign, bootstrap){

  stopifnot(inherits(psdesign, "psdesign"))
  stopifnot(inherits(bootstrap, "bootstrap"))

  psdesign <- bootstrap(psdesign)
  psdesign

}

#' Modify a psdesign object by adding on new components.
#'
#' This operator allows you to add objects to a psdesign object, such as integration models and risk models
#'
#' @param p1 An object of class \link{psdesign}
#' @param p2 Another object to be added to \code{p1}, see list below for possible options
#'
#' If the first object is an object of class \code{psdesign}, you can add
#' the following types of objects, and it will return a modified psdesign
#' object. Users will generally add them in the order that they appear.
#'
#' \itemize{
#'   \item \code{integration}: Add or replace integration model
#'   \item \code{riskmodel}: Add or replace risk model
#'   \item \code{estimate}: Estimate parameters
#'   \item \code{bootstrap}: Bootstrap estimates
#' }
#'
#' @export

"+.ps" <- function(p1, p2){

  if(inherits(p2, "integration")){
    add_integration(p1, p2)
    } else if(inherits(p2, "riskmodel")){
      add_riskmodel(p1, p2)
    } else if(inherits(p2, "estimate")){
      add_estimate(p1, p2)
    } else if(inherits(p2, "bootstrap")){
      add_bootstrap(p1, p2)
    }

}
