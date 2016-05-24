#' Manning's n for Corrugated Metal Pipe
#'
#' A data.frame containing the type of channel and description along with the
#' minimum, normal, and maximum value of n, if it exists. n is the "Gauckler-
#' Manning coefficient (commonly called Manning's n)" and it's dimensionless.
#' Source: Manning formula.
#'
#' @format A data frame with 25 rows and 2 variables:
#' \describe{
#' \item{Type of Pipe and Diameter and Corrugation Dimension}{Name of the type
#'    of conduit and any descriptive information}
#' \item{n}{Manning's n}
#' }
#'
#'
#' @references
#' \enumerate{
#'    \item This data is from FishXing Version 3.0 Beta (2006) by Michael Furniss, Michael Love, Susan Firor, Kathleen Moynan, Antonio Llanos, Jeff Guntle, and Robert Gubernick. See \url{http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm}. The original data source is Ven Te Chow, \emph{Open-Channel Hydraulics}, New York City, York: McGraw-Hill, 1959.
#'    \item Wikimedia Foundation, Inc. Wikipedia, 26 November 2015, “Manning formula”, \url{https://en.wikipedia.org/wiki/Manning_formula}.
#' }
#'
#'
#' @docType data
#' @name nmetalpipe
#' @usage nmetalpipe
#' @examples
#' nmetalpipe
NULL
