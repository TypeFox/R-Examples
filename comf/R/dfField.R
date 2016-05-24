#' Field data example
#'
#' Randomly sampled data from a field study campaign
#' with data from 156 samples. For further description, 
#' see the reference given.
#'
#' @docType data
#'
#' @usage data(dfField)
#'
#' @format A data frame with 156 rows and 9 variables:
#' \describe{
#'   \item{ta}{air temperature [degree C]}
#'   \item{tr}{radiant temperature [degree C] - same as ta}
#'   \item{rh}{relative humidity [\%]}
#'   \item{trrm}{running mean outdoor temperature [degree C]}
#'   \item{clo}{clothing insulation level [CLO]}
#'   \item{tout}{outdoor air temperature [degree C]}
#'   \item{av}{indoor air velocity [m/s]}
#'   \item{met}{metabolic rate [MET]}
#'   \item{asv}{actual thermal sensation vote on ASHRAE scale [ ]}
#' }
#'
#' @keywords datasets
#'
#' @references Schweiker, M. and Wagner, A. Exploring potentials and limitations of the adaptive thermal heat balance framework. Proceedings of 9th Windsor Conference: Making Comfort Relevant Cumberland Lodge, Windsor, UK, 2016.
#' (\href{http://www.windsorconference.com}{nceub})
#'
#' @examples
#' data(dfField)
#' head(dfField)
"dfField"