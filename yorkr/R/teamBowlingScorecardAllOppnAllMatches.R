##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBowlingScorecardAllOppnAllMatches
# This function computes the performance of bowlers of team against all opposition in all matches
# This function returns a dataframe
#
###########################################################################################
#' @title
#' Team bowling scorecard all opposition all matches
#'
#' @description
#' This function computes returns the bowling dataframe of bowlers deliveries, maidens, overs, wickets
#' against all oppositions in all matches
#'
#' @usage
#' teamBowlingScorecardAllOppnAllMatches(matches,theTeam)
#'
#' @param matches
#' The matches of the team against all oppositions and all matches
#'
#' @param theTeam
#' Team for which bowling performance is required
#'
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
#' matches <-getAllMatchesAllOpposition("India",dir="../data/",save=TRUE)
#'
#' # Or load directly from saved file
#' # load("allMatchesAllOpposition-India.RData")
#'
#' # Top opposition bowlers performances against India
#' teamBowlingScorecardAllOppnAllMatches(matches,"India")
#'
#' #Top Indian bowlers against respective opposition
#' teamBowlingScorecardAllOppnAllMatches(matches,'Australia')
#' teamBowlingScorecardAllOppnAllMatches(matches,'South Africa')
#' teamBowlingScorecardAllOppnAllMatches(matches,'England')
#' }
#'
#' @seealso
#' \code{\link{teamBowlingScorecardAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesPlot}}\cr
#'
#' @export
#'

teamBowlingScorecardAllOppnAllMatches <- function(matches,theTeam){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=NULL
    over=wickets=maidens=NULL
    a <-filter(matches,team==theTeam)

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
    i <- summarise(group_by(h,bowler),wickets=length(wicketPlayerOut))

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
