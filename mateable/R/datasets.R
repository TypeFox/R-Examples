#' Information about mating scene at site eelr in 2012.
#'
#' This dataframe contains information about all 44 plants that flowered in 2012
#' at the site eelr (East Elk Lake Road). Kelly Kapsar visited plants regularly
#' to determine the starting and ending dates of flowering on every head of
#' every plant. The metadata for the phenology dataset can be provide upon request
#' to the maintainer. Plants were mapped with gps with
#' better than 6 cm precision.
#'
#' @section Variables:
#' Variables:
#' \itemize{
#' \item tagNo, unique identifier for each plant
#' \item heads, number of flowering heads per plant in 2012
#' \item firstDay, the first day that any head on the plant shed pollen
#' \item lastDay, the last day that any head on the plant shed pollen
#' \item Ecoord, the x-coordinate of each plant in meters
#' \item Ncoord, the y-coordinate of each plant in meters
#' }
#' @docType data
#' @name eelr2012
#' @usage eelr2012
#' @format A 44 x 6 data frame
#' @references Wagenius, S. 2006. Scale-dependence of reproductive failure in
#' fragmented \emph{Echinacea} populations. Ecology 87: 931-941.
#' @keywords datasets
#' @examples
#' dim(eelr2012)
#' str(eelr2012)
NULL
