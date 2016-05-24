##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 22 Mar 2016
# Function: teamBatsmenPartnershiOppnAllMatches
# This function computes the batting partnership of a team in all matches against
# an opposition. The report generated can be detailed or a summary
#
###########################################################################################
#' @title
#' Team batting partnership against a opposition all matches
#'
#' @description
#' This function computes the performance of batsmen against all bowlers of an oppositions in all matches. This
#' function returns a dataframe
#'
#' @usage
#' teamBatsmenPartnershiOppnAllMatches(matches,theTeam,report="summary")
#'
#' @param matches
#' All the matches of the team against the oppositions
#'
#' @param theTeam
#' The team for which the the batting partnerships are sought
#'
#' @param report
#' If the report="summary" then the list of top batsmen with the highest partnerships is displayed. If
#' report="detailed" then the detailed break up of partnership is returned as a dataframe
#'
#' @return partnerships
#' The data frame of the partnerships
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
#' # Get all matches for team India against all oppositions
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#' # You can also directly load the data
#' #load("India-Australia-allMatches.RData")
#'
#' m <-teamBatsmenPartnershiOppnAllMatches(a,'India',report="summary")
#' m <-teamBatsmenPartnershiOppnAllMatches(a,'Australia',report="detailed")
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBowlersVsBatsmenMatch}}\cr
#' \code{\link{teamBattingScorecardMatch}}\cr
#'
#' @export
#'
teamBatsmenPartnershiOppnAllMatches <- function(matches,theTeam,report="summary"){

    team=batsman=nonStriker=partnershipRuns=runs=totalRuns=NULL
    a <-filter(matches,team==theTeam)
    #Get partnerships
    df <- data.frame(summarise(group_by(a,batsman,nonStriker),sum(runs)))
    names(df) <- c("batsman","nonStriker","partnershipRuns")
    b <- summarise(group_by(df,batsman),totalRuns=sum(partnershipRuns))
    c <- arrange(b,desc(totalRuns))
    d <- full_join(df,c,by="batsman")
    if(report == "detailed"){
        partnerships <- arrange(d,desc(totalRuns))
    } else{
        partnerships <- arrange(c,desc(totalRuns))
    }
    partnerships
}
