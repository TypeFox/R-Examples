#' Bayesian Approaches for HTE Analysis
#'
#' @docType package
#' @name beanz-package
#' @aliases beanz
#' @useDynLib beanz, .registration = TRUE
#'
#' @importFrom rstan sampling extract stanc
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines par plot points text
#'
#' @import survival
#' @import stats
#' @import Rcpp
#' @import methods
#'
#'
#' @description
#'
#' This package contains the functions for running Bayesian models implemented in
#' \code{STAN} for HTE analysis.
#'
#'
#' @section Notation:
#'
#' Consider a randomized two-arm clinical trial. Let \eqn{Y} denote the response
#'     and \eqn{Z} denote treatment arm assignment. For subgroup analysis,
#'     assume there are \eqn{P} baseline covariates, \eqn{X_1,\ldots,X_P}, of
#'     interest. The covariates can be binary, ordinal with numerical values, or
#'     nominal variables. Let \eqn{\Omega = \{(X_1,\ldots,X_P)\}} denote the
#'     collection of subgroups defined by the covariates. Let \eqn{\theta_g}
#'     denote the treatment effect in subgroup \eqn{G=g}, and let
#'     \eqn{\widehat{\theta}_g} be the estimated \eqn{\theta} in subgroup
#'     \eqn{G=g} with \eqn{\widehat{\sigma}^2_g} the estimated variance
#'     associated with \eqn{\widehat{\theta}_g}.
#'
#' @section Models:
#'
#' We approximate the distribution of \eqn{\widehat{\theta}_g}  by
#' \deqn{
#' \widehat{\theta}_g | \theta_g, \sigma^2_g \sim N(\theta_g, \sigma^2_g)
#' }
#' and assign an informative prior to \eqn{\sigma^2_g}
#'
#' \deqn{ \log \sigma_g^2 | \widehat{\sigma}^2_g, \Delta_1, \Delta_2 \sim
#' Unif( \log \widehat{\sigma}^2_g - \Delta_1, \log\widehat{\sigma}^2_g +
#' \Delta_2) }
#'
#' where \eqn{\Delta_1} and \eqn{\Delta_2} are parameters specified
#' by the users.
#'
#' We consider a set of models together with the priors for \eqn{\theta_g}:
#'
#' \describe{
#'
#' \item{\strong{No subgroup effect model}}{ This model assumes that
#'  patients in all the subgroups are exchangeable. That is, all the subgroups
#'  are statistically identical with regard to the treatment effect and there is
#'  no subgroup effect. Information about treatment effects can be directly
#'  combined from all subgroups for inference. The model is specified as follows:
#'
#' \deqn{
#' \begin{array}{rcl}
#'     \theta_g &=& \mu\\
#'      \mu     &\sim& N(0, B),
#' \end{array}
#' }
#'
#' where \eqn{B} is large in relation to the magnitude of the
#'  treatment effect size so that the prior for \eqn{\mu} is essentially
#'  non-informative. }
#'
#' \item{\strong{Full stratification model}}{ The subgroups are fully
#' distinguished from each other with regard to the treatment effect. There is
#' no information about treatment effects shared between any subgroups. The
#' model is specified as follows:
#'
#'\deqn{
#' \begin{array}{rcl}
#'   \theta_g &=& \mu_g \\
#' \mu_g    &\sim& N(0, B).
#' \end{array}
#' }
#'}
#'
#' \item{\strong{Simple regression model}}{ The model introduces a first-order,
#'  linear regression structure. This model takes into account the information
#'  that the subgroups are formulated based on the set of baseline covariates.
#'  The coefficients are assumed to be exchangeable among subgroups. Information
#'  about treatment effects are shared between subgroups with similar baseline
#'  covariates through these coefficients. The model is specified as follows:
#'
#'  \deqn{ \begin{array}{rcl}
#'\theta_g|X_g &=& \mu + \sum_{j=1}^P X'_{g,j} \gamma_j \\
#' \mu &\sim& N(0,B) \\
#' \gamma_j &\sim& N(0, B) \qquad j=1,\ldots,P.
#' \end{array} }
#' }
#'
#' \item{\strong{Basic shrinkage model}}{ This approach assumes all subgroups
#' are exchangeable with regards to the treatment effect. The model is specified
#' as follows:
#'
#' \deqn{
#' \begin{array}{rcl}
#' \theta_g  &=& \mu + \phi_g \\
#' \mu      &\sim& N(0, B) \\
#' \phi_g   &\sim& N(0, \omega^2) \\
#' \omega &\sim& {Half-}N(C).
#' \end{array} }
#' }
#'
#'
#' \item{\strong{Simple regression and shrinkage model}}{
#' This model combines basic regression with shrinkage, with a linear
#' regression structure and a random effect term. Direct estimates are shrunken
#' towards the regression surface. The model is specified as follows:
#'
#' \deqn{ \begin{array}{rcl}
#' \theta_g   &=& \mu + \sum_{j=1}^P  X'_{g,j} \gamma_j + \phi_g \\
#' \mu &\sim& N(0,B) \\
#' \gamma_j &\sim& N(0, 1 B) \qquad j=1,\ldots,P\\
#' \phi_g   &\sim& N(0, \omega^2) \\
#' \omega &\sim& {Half-}N(C).
#'
#' \end{array}
#' } }
#'
#'
#' \item{\strong{Dixon and Simon model}}{ This model assumes that the elements in coefficient
#' are exchangeable with each other, which allows information sharing among
#' covariate effects. Similar to the simple regression model, only the first-order
#' interactions are considered. The model is specified as follows:
#'
#' \deqn{
#' \begin{array}{rcl}
#'
#' \theta_g   &=& \mu + \sum_{j=1}^P X'_{g,j} \gamma_j \\
#' \mu        &\sim& N(0,B) \\
#' \gamma_j   &\sim& N(0, \omega^2) \\
#' \omega   &\sim& {Half-}N(C).
#'
#' \end{array} } }
#'
#' \item{\strong{Extended Dixon and Simon model}}{
#' This approach extends the Dixon and Simon model by introducing the
#' higher-order interactions, with the interaction effects exchangeable.
#' The model is specified as follows:
#'\deqn{
#' \begin{array}{rcl}
#'  \theta_g   &=& \mu + \sum_{k=1}^P \sum_{j \in \xi^{(k)}}  X'_{\xi^{(k)},j} \gamma^{(k)}_{j} \\
#' \mu &\sim& N(0,B) \\
#' \gamma^{(k)}_j &\sim& N(0, \omega_k^2) \qquad k=1,\ldots,P, \quad j\in \xi^{(k)} \\
#' \omega_k &\sim& {Half-}N(C),
#' \end{array}
#' }
#' where \eqn{\xi^{(k)}} denotes the set of \eqn{k}th order interaction terms
#'}
#'
#'
#'}
#'
#'@section Graphical user interface (GUI):
#'
#' This package provides a web-based Shiny GUI. See \code{\link{run.beanz}} for
#' details.
#'
#' @references
#'
#' Jones HE, Ohlssen DI, Neuenschwander B, Racine A, Branson M (2011). Bayesian
#' models for subgroup analysis in clinical trials. Clinical Trials, 8(2),
#' 129-143.
#'
#' Dixon DO, Simon R (1991). Bayesian subset analysis. Biometrics, 47(3), 871-881.
#'
NULL
