#' Cervical cancer data
#' 
#' The data set contains the number of cervical cancer deaths (ICD 180) and 
#' woman-years at risk for four age groups in four different European countries 
#' during 1969-1973.
#' 
#' @docType data
#' @usage data(cervical)
#' @format A data frame with 16 rows and 19 variables: 
#' \describe{
#'  \item{\code{y}}{number of cervical cancer deaths for different age categories
#'    and European countries between 1969-1973}
#'  \item{\code{E}}{number of woman-years at risk (given in thousands)}
#'  \item{\code{country}}{factor variable of European countries}
#'  \item{\code{agegroup}}{factor variable of age categories}
#'  \item{\code{X.1}, \code{X.2}, \code{X.3}}{predictor variables for country effects
#'    using dummy coding (i.e. England, France, Italy)}
#'  \item{\code{X.4}, \code{X.5}, \code{X.6}}{predictor variables for age effects
#'    using dummy coding (i.e. 35-44, 45-54, 55-64, in years)}
#'  \item{\code{X.7}, \code{X.8}, \code{X.9}, \code{X.10}, \code{X.11}, \code{X.12}, \code{X.13}, \code{X.14}, \code{X.15}}{predictor
#'    variables for interaction effects of age and country}
#' }
#' 
#' @note The lowest age category (25-34) in Belgium is used as the reference category.
#' @source World Health Organization (1976). World Health Statistics Annual: 1969-1976,
#'  Vol. I, \emph{Vital Statistics and Causes of Death}. Geneva: WHO.
#' @source Whittemore, A. S. and Gong, G. (1991). Poisson regression 
#'  with missclassified counts: Application to cervical cancer mortality rates. 
#'  \emph{Applied Statistics}, \strong{40}, 81-93. 
#' @seealso \code{\link{cervical_validation}}, \code{\link{pogitBvs}}
#' @name cervical
#' @keywords datasets
NULL


