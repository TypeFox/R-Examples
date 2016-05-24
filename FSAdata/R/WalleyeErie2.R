#' @title Biological data for Walleye from Lake Erie, 2003-2014.
#' 
#' @description Walleye (\emph{Sander vitreus}) biological data (length, weight, sex, maturity, and age) from several locations in Lake Erie, October-November of 2003-2014.
#' 
#' @name WalleyeErie2
#' @docType data
#' @format A data frame of 33734 observations on the following 10 variables:
#' \describe{
#'  \item{setID}{Unique gear (multifilament gill net kegged 6 ft below surface) set identification number.} 
#'  \item{loc}{Regional location (1=Toledo to Huron, 2=Huron to Fairport Harbor, 3=Fairport Harbor Conneaut).} 
#'  \item{grid}{2.5-minute sampling grid location}
#'  \item{year}{Year of data}
#'  \item{tl}{Total length (mm)}
#'  \item{w}{Weight (g).  There are several missing values}
#'  \item{sex}{Sex (female, male)}
#'  \item{mat}{Maturity (immature, mature)}
#'  \item{age}{Age (yrs) from otoliths}
#' }
#' @section Topic(s): \itemize{
#'  \item Growth
#'  \item von Bertalanffy
#'  \item Weight-Length
#'  \item Catch curve
#'  \item Mortality
#'  \item Maturity
#'  \item Size Structure
#'  \item Length Frequency
#'  \item Condition
#' }
#' @concept 'Growth' 'Weight-Length' 'Catch curve' 'Mortality' 'Maturity' 'Size Structure' 'Condition' 'Length Frequency'
#' @source These unpublished data are from the Ohio Department of Natural Resources, Division of Wildlife (via Christopher Vandergoot).  \bold{Do not use for other than educational purposes without permission from the source.}
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeErie2)
#' str(WalleyeErie2)
#' head(WalleyeErie2)
#' xtabs(~year+loc+sex,data=WalleyeErie2)
#' 
NULL