#' MF Package
#'
#' Includes functions related to mitigated fraction. \cr \cr For internal use only at
#' the USDA Center for Veterinary Biologics. \cr
#'
#' \tabular{ll}{
#' Package: \tab MF-package\cr
#' Type: \tab Package\cr
#' Version: \tab 4.3.2\cr
#' Date: \tab 2014-01-10\cr
#' License: \tab MIT \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name MF-package
#' @aliases MF
#' @docType package
#' @author David Siev \email{David.Siev@@aphis.usda.gov}
#' @examples
#' #---------------------------------------------
#' # Checking MF package
#' #---------------------------------------------
#' example(MFr)
#' #---------------------------------------------
#' # End examples
#' #---------------------------------------------
#' invisible()
NA



#' @name calflung
#' @title calflung dataset
#' @alias calflung-data
#' @docType data
#' @description something here
#' @format a data frame with 50 observations of the following 2 variables, no NAs
#' \describe{
#' \item{group}{Treatment group. One of con = control or vac = vaccinate }
#' \item{lesion}{Percent lung lesion, in decimal form}
#' }
#' @keywords datasets
NA


#' @name mlesions
#' @title mlesions dataset
#' @alias mlesions-data
#' @docType data
#' @description something here
#' @format a data frame with 52 observations of the following 3 variables, no NAs
#' \describe{
#' \item{cage}{Cage ID. 1 - 26}
#' \item{tx}{Treatment. One of 'con' or 'vac'}
#' \item{les}{Percent lung lesion}
#' }
#' @keywords datasets
NA

#' @name piglung
#' @title piglung dataset
#' @alias piglung-data
#' @docType data
#' @description something here
#' @format a data frame with 102 observations of the following 3 variables, no NAs
#' \describe{
#' \item{lesion}{Percent lung lesion}
#' \item{group}{Treatment group. One of 'con' or 'vac'}
#' \item{litter}{Litter ID}
#' }
#' @keywords datasets
NA