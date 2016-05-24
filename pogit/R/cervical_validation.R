#' Cervical cancer valiation data
#' 
#' Additionally to the main study sample (see \code{\link{cervical}}), 
#' validation data are available that give information on how likely physicians 
#' from different countries are to identify and correctly report a true cervical 
#' cancer death. In that study, a sample of physicians in each country completed 
#' a death certificate for one specific patient who had died of cervical cancer 
#' and the number of correct death certificates in each country was recorded. 
#' Validation data are therefore available on country level but provide no 
#' information on the reporting probability specific for age. 
#' 
#' @docType data
#' @usage data(cervical_validation)
#' @format A data frame with 4 rows and 6 variables: 
#' \describe{ 
#'  \item{\code{v}}{number of correct death certificates in each country 
#'    in the validation sample } 
#'  \item{\code{m}}{size of validation sample in each country} 
#'  \item{\code{country}}{factor variable of European countries} 
#'  \item{\code{W.1}, \code{W.2}, \code{W.3}}{predictor variables for country effects 
#'    using dummy coding (i.e. England, France, Italy)} 
#'  }
#'   
#' @note Belgium is used as the reference category.  
#' @source Kelson, M. and Farebrother, M. (1987). The Effect of Inaccuracies in 
#'   Death Certification and Coding Practices in the European Economic Community
#'   (EEC) on International Cancer Mortality Statistics. \emph{International
#'   Journal of Epidemiology}, 16, \strong{3}, 411-414.
#' @source Whittemore, A. S. and Gong, G. (1991). Poisson regression with
#'   missclassified counts: Application to cervical cancer mortality rates. 
#'   \emph{Applied Statistics}, \strong{40}, 81-93.
#' @seealso \code{\link{cervical}}, \code{\link{pogitBvs}}
#' @name cervical_validation
#' @keywords datasets
NULL

