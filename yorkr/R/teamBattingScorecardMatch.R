##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 20 Mar 2016
# Function: teamBattingScorecardMatch
# This function gets the batting scorecard of team in a match. The result is
# returned as a data frame
#
###########################################################################################
#' @title
#' Team batting scorecard of a team in a match
#'
#' @description
#' This function computes returns the batting scorecard (runs, fours, sixes, balls played) for the
#' team
#' @usage
#' teamBattingScorecardMatch(match,theTeam)
#'
#' @param match
#' The match for which the score card is required e.g.
#'
#' @param theTeam
#' Team for which scorecard required
#'
#' @return scorecard
#' A data frame with the batting scorecard
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
#' @examples
#' \dontrun{
#' a <- getMatchDetails("England","Pakistan","2006-09-05",dir="../temp")
#' teamBowlingScorecardMatch(a,'England')
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#'
#' @export
#'

teamBattingScorecardMatch <- function(match,theTeam){
    team=batsman=runs=fours=sixes=NULL
    byes=legbyes=noballs=wides=NULL
    a <-filter(match,team==theTeam)
    sz <- dim(a)
    if(sz[1] == 0){
        cat("No batting records.\n")
        return(NULL)
    }
    b <- select(a,batsman,runs)

    names(b) <-c("batsman","runs")
    #Compute the number of 4s
    c <-
        b %>%
        mutate(fours=(runs>=4 & runs <6)) %>%
        filter(fours==TRUE)

    # Group by batsman. Count 4s
    d <-    summarise(group_by(c, batsman),fours=n())

    # Get the total runs for each batsman
    e <-summarise(group_by(a,batsman),sum(runs))
    names(b) <-c("batsman","runs")
    details <- full_join(e,d,by="batsman")
    names(details) <-c("batsman","runs","fours")

    f <-
        b %>%
        mutate(sixes=(runs ==6)) %>%
        filter(sixes == TRUE)

    # Group by batsman. oOunt 6s
    g <- summarise(group_by(f, batsman),sixes=n())
    names(g) <-c("batsman","sixes")
    #Full join with 4s and 6s
    details <- full_join(details,g,by="batsman")

    # Count the balls played by the batsman
    ballsPlayed <-
        a  %>%
        select(batsman,byes,legbyes,wides,noballs,runs) %>%
        filter(wides ==0,noballs ==0,byes ==0,legbyes == 0) %>%
        select(batsman,runs)
    ballsPlayed<- summarise(group_by(ballsPlayed,batsman),count=n())
    names(ballsPlayed) <- c("batsman","ballsPlayed")


    details <- full_join(details,ballsPlayed,by="batsman")
    cat("Total=",sum(details$runs),"\n")
    # If there are NAs then
    if(sum(is.na(details$fours)) != 0){
         details[is.na(details$fours),]$fours <- 0
    }
    if(sum(is.na(details$sixes)) != 0){
        details[is.na(details$sixes),]$sixes <- 0
    }
    # Out the details
    details <- select(details,batsman,ballsPlayed,fours,sixes,runs)
    details

}
