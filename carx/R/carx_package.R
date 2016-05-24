#' \code{carx}: A package to fit Censored Auto-Regressive model with eXogenous covariates (CARX)
#'
#' \code{carx} is a package to estimate the parameters of the Censored AutoRegressive model with eXogenous
#' covariates (CARX), which can also be viewed as regression models with
#' censored responses
#' and autoregressive residuals. \code{carx} allows left, right, or interval
#' censoring for the response variable. The regression errors are assumed to
#' follow an autoregressive model with normal innovations. 
#' In addition to the estimation method, the package also
#' contains functions to predict future values, diagnose whether the model is
#' adequate, and plot functions to illustrate the data and model.
#'
#' More specifically, we estimate the parameters assumed in the following model.
#' Let \eqn{(Y_t)} be a censored time series with the latent process denoted by
#' \eqn{(Y_t^*)}.
#' For each \eqn{Y_t^*}, it can be censored by
#' either \eqn{(-\infty,c_{l,t})} or
#' \eqn{(c_{u,t},\infty)}, and if it is censored, \eqn{Y_t} will be recorded as \eqn{c_{l,t}}
#' or \eqn{c_{u,t}} respectively.
#'
#' The latent process \eqn{(Y_t^*)} is modelled as
#' \deqn{
#' Y_t^* = X_t' \beta + \eta_t,
#' }
#' and
#' \deqn{
#' \eta_t = \sum_{i=1}^p \psi_i \eta_{t-i} + \varepsilon_t,
#' }
#' where \eqn{(X_t)} is a covariate process with all values observable,
#' and the innovations \eqn{(\varepsilon_t)} are independent and identically normally distributed with mean 0 and variance \eqn{\sigma^2}.
#'
#' In this package we implemented the quasi-maximum likelihood estimator proposed by Wang and Chan (2015), for more details, please refer to the paper.
#'
#'
#'
#' @references Wang C, Chan KS (2015). "Quasi-likelihood estimation of a censored autoregressive model with exogenous variables." Submitted.
#'
#' @docType package
#' @name carx
#' @importFrom stats AIC coef fitted formula window tsdiag residuals
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL
