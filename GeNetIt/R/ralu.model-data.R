#' @name ralu.model
#' @docType data
#'
#' @title Columbia spotted frog (Rana luteiventris) data for specifying gravity model. Note, the data.frame is already log transformed. 
#'
#' @description Subset of data used in Murphy et al., (2010)
#'
#' @format A data.frame with 190 rows (sites) and 19 columns (covariates):
#' \describe{
#' \item{ARMI_ID}{Unique ID}
#' \item{FROM_SITE}{Unique from site ID}
#' \item{TO_SITE}{Unique to site ID}
#' \item{FST}{FST genetic distance}
#' \item{DPS}{DPS genetic distance}
#' \item{DISTANCE}{Graph edge distance}
#' \item{DEPTH_F}{At site water depth}
#' \item{HLI_F}{Heat Load Index}
#' \item{CTI_F}{Wetness Index}
#' \item{DEPTH_T}{At site water depth}
#' \item{HLI_T}{Heat Load Index}
#' \item{CTI_T}{Wetness Index}
#' \item{hli}{Heat Load Index}
#' \item{cti}{Wetness Index}
#' \item{ffp}{Frost Free Period}
#' \item{err27}{Roughness at 27x27 scale}
#' \item{rsp}{Relative Slope Position}
#' \item{ridge}{Percent Ridge Line}
#' \item{hab_ratio}{Ratio of suitable dispersal habitat}
#' }
#'
#' @references
#' Murphy M.A., R. Dezzani, D.S. Pilliod & A.S. Storfer (2010) Landscape genetics of high mountain frog metapopulations. Molecular Ecology 19(17):3634-3649 
#'
NULL
