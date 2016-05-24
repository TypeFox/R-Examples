#' @docType data
#'
#' @name salmonellaSerologyData
#'
#' @title
#' Salmonella Antibody Levels Data
#'
#' @description
#' Antibody levels measured in cross-sectional population sample for Salmonella with SSI Mixed ELISA procedure.
#'
#' @usage
#' data("salmonellaSerologyData")
#'
#' @format
#' A data frame with 456 observations on the following 6 variables:
#' \describe{
#' \item{\code{id}}{Observation id}
#' \item{\code{IgG}}{Measured IgG antibody level}
#' \item{\code{IgM}}{Measured IgM antibody level}
#' \item{\code{IgA}}{Measured IgA antibody level}
#' \item{\code{sex}}{Sex indicator (1 - Male, 2 - Female)}
#' \item{\code{age}}{Age}
#' }
#'
#' @source
#' Dutch Pienter 1 serosurvey
#'
#' @references
#' Teunis, P. F., van Eijkeren, J. C., Ang, C. W., van Duynhoven, Y. T., Simonsen, J. B., Strid, M. A., van Pelt W.\cr
#' "Biomarker Dynamics: Estimating Infection Rates From Serological Data"\cr
#' Statistics in Medicine 31, no. 20 (September 9, 2012): 2240--48. doi:10.1002/sim.5322.
#'
#' @examples
#'
#' # show first rows of the data
#' head(salmonellaSerologyData)
#'
#' # plot the data
#' plot(salmonellaSerologyData)
#'
NULL
