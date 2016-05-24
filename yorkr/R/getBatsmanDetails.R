##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: getBatsmanDetails
# This function gets the batting details of a batsman
#
###########################################################################################
#' @title
#' Get batting details of  batsman from match
#'
#' @description
#' This function gets the batting details of a batsman given the match data as a RData file
#' @usage
#' getBatsmanDetails(team,name,dir=".")
#'
#' @param team
#' The team of the batsman e.g. India
#'
#' @param name
#' Name of batsman
#'
#' @param dir
#' The directory where the source file exists
#'
#' @return None
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#'
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#'
#' @examples
#' \dontrun{
#' getBatsmanDetails(team="India",name="Kohli",dir=pathToFile)
#' }
#'
#'
#' @seealso
#' \code{\link{batsmanRunsPredict}}\cr
#' \code{\link{batsmanMovingAverage}}\cr
#' \code{\link{bowlerWicketsVenue}}\cr
#' \code{\link{bowlerMeanRunsConceded}}\cr
#'
#' @export
#'
#'

getBatsmanDetails <- function(team, name,dir="."){
    batsman=battingDetails=NULL
    fl <- paste(dir,"/",team,"-BattingDetails.RData",sep="")
    print(fl)
    load(fl)
    details <- battingDetails
    batsmanDetails <- filter(details,grepl(name,batsman))
    batsmanDetails
}
