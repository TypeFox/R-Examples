##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: bowlerDeliveryWickets
# This function creates a data frame of deliveries bowled and wickets
#
###########################################################################################
#' @title
#' Number of deliveries to wickets
#'
#' @description
#' This function creates a dataframe of balls bowled versus the wickets taken by
#' the bowler
#' @usage
#' bowlerDeliveryWickets(match,theTeam,name)
#'
#' @param match
#' Data frame of the match
#'
#' @param theTeam
#' The team for which the delivery wickets have to be computed
#'
#' @param name
#' The name of the bowler
#'
#' @return dataframe
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' #Get match data
#' match <- getMatchDetails("England","Pakistan","2006-09-05",dir="../data")
#' bowlerDeliveryWickets(match,"India","Jadeja")
#' }
#'
#' @seealso
#' \code{\link{batsmanFoursSixes}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#' \code{\link{bowlerDeliveryWickets}}\cr
#' \code{\link{bowlerMeanEconomyRate}}\cr
#' \code{\link{bowlerMeanRunsConceded}}\cr
#'
#'

bowlerDeliveryWickets <- function(match,theTeam,name){
    team = bowler = wicketPlayerOut =delivery = wicketNo = NULL
    d <- NULL
    a <-filter(match,team!=theTeam)
    b <- filter(a,grepl(name,bowler))
    if(dim(b)[1] != 0){
        b$delivery<- seq(1:dim(b)[1])
        c <- filter(b,wicketPlayerOut != "nobody")
        if(dim(c)[1] !=0){
            c$wicketNo <- seq(1:dim(c)[1])
            d <- select(c,bowler,delivery,wicketNo,date)
        }

    }
    d
}
