##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 21 Mar 2016
# Function: teamBowlingScorecardMatch
# This function computes the performance of bowlers of team in a match.
# This function returns a dataframe
#
###########################################################################################
#' @title
#' Compute and return the bowling scorecard of a team in a match
#'
#' @description
#' This function computes and returns the bowling scorecard of a team in a match
#'
#' @usage
#' teamBowlingScorecardMatch(match,theTeam)
#'
#' @param match
#' The match between the teams
#'
#' @param theTeam
#' Team for which bowling performance is required
#'
#' @return l
#' A data frame with the bowling performance in alll matches against all oppositions
#'
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' #Get all matches between India  and other opposition
#' a <- getMatchDetails("England","Pakistan","2006-09-05",dir="../temp")
#'
#' teamBowlingScorecardMatch(a,'England')
#' }
#'
#' @seealso
#' \code{\link{teamBowlingWicketMatch}}\cr
#' \code{\link{teamBowlersVsBatsmenMatch}}\cr
#' \code{\link{teamBattingScorecardMatch}}\cr
#'
#' @export
#'
teamBowlingScorecardMatch <- function(match,theTeam){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=over=NULL
    # Compute the maidens,runs conceded and overs for the bowlers.
    # The bowlers performance of the team is got when the other side is batting. Hence '!-"
    a <-filter(match,team != theTeam)
    a1 <- unlist(strsplit(a$ball[1],"\\."))
    # Create a string for substitution 1st or 2nd
    a2 <- paste(a1[1],"\\.",sep="")
    # only wides and noballs need to be included with runs for bowlers.
    # Note: byes and legbyes should not be included
    b <-  a %>%
        select(bowler,ball,noballs,wides,runs,wicketKind,wicketPlayerOut) %>%
        #mutate(over=gsub("1st\\.","",ball)) %>%
        mutate(over=gsub(a2,"",ball)) %>%
        mutate(over=gsub("\\.\\d+","",over))

    #Calculate the number of maiden overs
    c <- summarise(group_by(b,bowler,over),sum(runs,wides,noballs))
    names(c) <- c("bowler","over","runsConceded")
    d <-summarize(group_by(c,bowler),maidens=sum(runsConceded==0))

    #Compute total runs conceded (runs_wides+noballs)
    e <- summarize(group_by(c,bowler),runs=sum(runsConceded))

    # Calculate the number of overs bowled by each bwler
    f <- select(c,bowler,over)
    g <- summarise(group_by(f,bowler),overs=length(unique(over)))


    #Compute number of wickets
    h <- b %>%
        select(bowler,wicketKind,wicketPlayerOut) %>%
        filter(wicketPlayerOut != "nobody")
    i <- summarise(group_by(h,bowler),wickets=length(unique(wicketPlayerOut)))

    #Join the over & maidens
    j <- full_join(g,d,by="bowler")
    # Add runs
    k <- full_join(j,e,by="bowler")
    # Add wickets
    l <- full_join(k,i,by="bowler")

    # Set NAs to 0
    if(sum(is.na(l$wickets)) != 0){
        l[is.na(l$wickets),]$wickets=0
    }
    l
}
