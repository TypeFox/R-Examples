#' @title
#' Simulate Cross-Sectional Serology Data
#'
#' @description
#' Function to simulate cross-sectional serology data for three antibodies: IgG, IgM, IgA.
#' The underlying parametric models are log-Normal distributions with parameters:
#' \itemize{
#' \item{IgG: \code{mu} = 1.246, \code{sd} = 0.757}
#' \item{IgM: \code{mu} = -0.764, \code{sd} = 0.413}
#' \item{IgA: \code{mu} = -1.13, \code{sd} = 0.318}
#' }
#'
#' @param
#' n Number of simulations, default 500.
#'
#' @return
#' A dataframe of simulated antibody data, with named columns IgG, IgM and IgA.
#'
#' @examples
#'
#' # simulate 500 observations of cross-sectional serology data
#' simulateSerologyData()
#'
#' # simulate 10 observations of cross-sectional serology data
#' simulateSerologyData(n = 10)
#'
#' @export
simulateSerologyData <- function(n = 500)
{
    simSerology <- function(n = 500, mu = 1, sd = 0.5)
    {
        return(exp(stats::rnorm(n = n, mean = mu, sd = sd)))
    }

    result <- data.frame(
        IgG = simSerology(n = n, mu = 1.246, sd = 0.757),
        IgM = simSerology(n = n, mu = -0.764, sd = 0.413),
        IgA = simSerology(n = n, mu = -1.13, sd = 0.318))

    return(result)
}
