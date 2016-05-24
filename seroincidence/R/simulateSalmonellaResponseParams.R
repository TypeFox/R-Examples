#' @title
#' Simulate Longitudinal Response Parameters for Salmonella
#'
#' @description
#' Simulate longitudinal response parameters \code{A} and \code{k} per antibody
#' for Salmonella (SSI Mixed ELISA).
#'
#' The underlying parametric models for peak levels \code{A} are Gamma distributions with parameters:
#' \itemize{
#' \item{IgG: \code{shape} = 1.1175, \code{scale} = 0.848}
#' \item{IgM: \code{shape} = 1.337, \code{scale} = 0.902}
#' \item{IgA: \code{shape} = 1.205, \code{scale} = 0.651}
#' }
#'
#' The underlying parametric models for decay rates \code{k} are Inverse Gamma distributions with parameters:
#' \itemize{
#' \item{IgG: \code{shape} = 0.869, \code{scale} = 1/728.5}
#' \item{IgM: \code{shape} = 0.731, \code{scale} = 1/514.6}
#' \item{IgA: \code{shape} = 1.759, \code{scale} = 1/132/5}
#' }
#'
#' @param
#' n Number of Monte Carlo simulations per antibody, default 500.
#'
#' @return
#' A list of two dataframes named \code{A} and \code{k}:
#' \describe{
#' \item{\code{A}}{Dataframe containing peak levels (SSI mixed ELISA units/ml). Named
#' columns contain estimated peak levels for IgG, IgM and IgA antibodies.}
#' \item{\code{k}}{Dataframe containing decay rates (1/days). Named columns contain
#' estimated decay rates for IgG, IgM and IgA antibodies.}
#' }
#'
#' @examples
#'
#' # simulate 500 observations of longitudinal response data
#' simulateSalmonellaResponseParams()
#'
#' # simulate 100 observations of longitudinal response data
#' simulateSalmonellaResponseParams(n = 100)
#'
#' @export
simulateSalmonellaResponseParams <- function(n = 500)
{
    result <- list(
        A = data.frame(
            IgG = stats::rgamma(n = n, shape = 1.175, scale = 0.848),
            IgM = stats::rgamma(n = n, shape = 1.337, scale = 0.902),
            IgA = stats::rgamma(n = n, shape = 1.205, scale = 0.651)),
        k = data.frame(
            IgG = .rinvgamma(n = n, shape = 0.869, scale = 1/728.5),
            IgM = .rinvgamma(n = n, shape = 0.731, scale = 1/514.6),
            IgA = .rinvgamma(n = n, shape = 1.759, scale = 1/132.5)))

    return(result)
}
