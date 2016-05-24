#' Improved Prediction Intervals for ARIMA Processes and Structural Time Series
#'
#' Package tsPI computes prediction intervals for ARIMA and structural
#' time series models by using importance sampling approach with uninformative priors
#' for model parameters, leading to more accurate coverage probabilities in frequentist sense.
#' Instead of sampling the future observations and hidden states
#' of the state space representation of the model, only model parameters are sampled,
#' and the method is based solving the equations corresponding to the conditional
#' coverage probability of the prediction intervals. This makes method relatively
#' fast compared to for example MCMC methods, and standard errors of prediction
#' limits can also be computed straightforwardly.
#'
#' @docType package
#' @name tsPI
#' @aliases tsPI
#' @useDynLib tsPI
#' @importFrom KFAS SSModel SSMtrend SSMarima SSMseasonal KFS simulateSSM
#' @importFrom stats arima arima.sim dnorm end frequency logLik optim pnorm predict rchisq rnorm sd toeplitz ts uniroot window na.exclude
#' @references
#' \enumerate{
#' \item{Helske, J. and Nyblom, J. (2013). Improved frequentist prediction
#' intervals for autoregressive models by simulation.
#' In Siem Jan Koopman and Neil Shephard, editors,
#' Unobserved Components and Time Series Econometrics. Oxford University Press. In press.}
#'  \item{Helske, J. and Nyblom, J. (2014). Improved frequentist prediction intervals for
#'  ARMA models by simulation.
#'  In Johan Knif and Bernd Pape, editors,
#' Contributions to Mathematics, Statistics, Econometrics, and Finance:
#' essays in honour of professor Seppo Pynnönen,
#' number 296 in Acta Wasaensia, pages 71–86. University of Vaasa.}
#' \item{Helske, J. (2015). Prediction and interpolation of time series by state space models. University of Jyväskylä. PhD thesis, Report 152. }
#' }
NULL
