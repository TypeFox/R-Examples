##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 23 Mar 2016
# Function: teamBattingBattingScorecardOppnAllMatches
# This function computes the batting scorecard for the team against an oppositon
# in all matches against the opposition
#
#
###########################################################################################
#' @title
#' Team batting scorecard of a team in all matches against an opposition
#'
#' @description
#' This function computes returns the batting scorecard (runs, fours, sixes, balls played) for the
#' team in all matches against an opposition
#'
#' @usage
#' teamBattingScorecardOppnAllMatches(matches,main,opposition)
#'
#' @param matches
#' the data frame of all matches between a team and an opposition obtained with
#' the call getAllMatchesBetweenteam()
#'
#' @param main
#' The main team for which scorecard required
#'
#' @param opposition
#' The opposition team
#'
#' @return scorecard
#' The scorecard of all the matches
#'
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
#' # Get all matches between India and Australia
#' matches <- getAllMatchesBetweenTeams("India","Australia",dir="../data",save=TRUE)
#' # Compute the scorecard of India in matches with australia
#' teamBattingScorecardOppnAllMatches(matches,main="India",opposition="Australia")
#'
#' #Get all matches between Australia and India
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#' #Compute the batting scorecard of Australia
#' teamBattingScorecardOppnAllMatches(matches,"Australia","India")
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#'
#' @export
#'
teamBattingScorecardOppnAllMatches <- function(matches,main,opposition){
    team=batsman=runs=fours=sixes=NULL
    byes=legbyes=noballs=wides=NULL

    a <-filter(matches,team==main)
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

    # Compute sixes
    f <-
        b %>%
        mutate(sixes=(runs ==6)) %>%
        filter(sixes == TRUE)
    # Group by batsman. COunt 6s
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
    details <- arrange(details,desc(runs),desc(sixes),desc(fours))
    details <- select(details,batsman,ballsPlayed,fours,sixes,runs)
    details

}
