##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 23 Mar 2016
# Function: teamBowlingPerfOppnAllMatches
# This function computes the team bowling performance against an opposition and picks
# the top bowlers based on number of wickets taken
#
#
###########################################################################################
#' @title
#' team bowling performance all matches against an opposition
#'
#' @description
#' This function computes returns the bowling dataframe of bowlers deliveries, maidens, overs, wickets
#' against an opposition in all matches
#'
#' @usage
#' teamBowlingPerfOppnAllMatches(matches,main,opposition)
#'
#' @param matches
#' The matches of the team against an opposition.
#'
#' @param main
#' Team for which bowling performance is required
#'
#' @param opposition
#' The opposition Team
#'
#'
#' @return l
#' A data frame with the bowling performance
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
#' #Get all matches between India and Autralia
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#'
#' # Or load directly from saved file
#' # load("India-Australia-allMatches.RData")
#'
#' teamBowlingPerfOppnAllMatches(matches,"India","Australia")
#' teamBowlingPerfOppnAllMatches(matches,main="Australia",opposition="India")
#' }
#'
#' @seealso
#' \code{\link{teamBowlersWicketsOppnAllMatches}}\cr
#' \code{\link{teamBowlersWicketRunsOppnAllMatches}}\cr
#' \code{\link{teamBowlersWicketKindOppnAllMatches}}\cr
#'
#' @export
#'

teamBowlingPerfOppnAllMatches <- function(matches,main,opposition){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=over=NULL
    wickets=maidens=NULL
    # Compute the maidens,runs conceded and overs for the bowlers
    a <-filter(matches,team != main)

    # only wides and noballs need to be included with runs for bowlers.
    # Note: byes and legbyes should not be included
    b <-  a %>%
        select(bowler,ball,noballs,wides,runs,wicketKind,wicketPlayerOut) %>%
        mutate(over=gsub("1st\\.","",ball)) %>%
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

    # Set NAs to 0 if there are any
    if(sum(is.na(l$wickets)) != 0){
        l[is.na(l$wickets),]$wickets=0
    }
    # Arrange in descending order of wickets and runs and ascending order for maidens
    l <-arrange(l,desc(wickets),desc(runs),maidens)
    l
}
