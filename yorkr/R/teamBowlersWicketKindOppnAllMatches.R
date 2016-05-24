##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 24 Mar 2016
# Function: teamBowlersWicketKindOppnAllMatches
# This function computes the the types of wickets taken by the bowlers, caught,bowled, c&b etc
#
#
###########################################################################################
#' @title
#' Team bowlers wicket kind against an opposition in all matches
#'
#' @description
#' This function computes performance of bowlers of a team and the wicket kind against an
#' opposition in all matches against the opposition
#'
#' @usage
#' teamBowlersWicketKindOppnAllMatches(matches,main,opposition,plot=TRUE)
#'
#' @param matches
#' The data frame of all matches between a team the opposition. This dataframe can be obtained with
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#'
#' @param main
#' The team for which the performance is required
#'
#' @param opposition
#' The opposing team
#'
#' @param plot
#' If plot=TRUE then a plot is displayed else a dataframe is returned
#'
#' @return None or dataframe
#' The return depends on the value of the plot
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
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#'
#' teamBowlersWicketKindOppnAllMatches(matches,"India","Australia",plot=TRUE)
#' m <- teamBowlersWicketKindOppnAllMatches(matches,"Australia","India",plot=FALSE)
#'
#' teamBowlersWicketKindOppnAllMatches(matches,"Australia","India",plot=TRUE)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesRept}}
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesPlot}}
#'
#' @export
#'
teamBowlersWicketKindOppnAllMatches <- function(matches,main,opposition,plot=TRUE){
    team=bowler=ball=NULL
    runs=over=runsConceded=NULL

    byes=legbyes=noballs=wides=runConceded=NULL
    extras=wicketFielder=wicketKind=wicketPlayerOut=NULL
    # Compute the maidens,runs conceded and overs for the bowlers
    a <-filter(matches,team !=main)

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


    #Compute number of wickets
    h <- b %>%
        select(bowler,wicketKind,wicketPlayerOut) %>%
        filter(wicketPlayerOut != "nobody")
    i <- summarise(group_by(h,bowler),wickets=length(unique(wicketPlayerOut)))

    r <- full_join(h,e,by="bowler")

    # Set NAs to 0
    if(sum(is.na(r$wicketKind)) != 0){
        r[is.na(r$wicketKind),]$wicketKind="noWicket"
    }
    if(sum(is.na(r$wicketPlayerOut)) !=0){
        r[is.na(r$wicketPlayerOut),]$wicketPlayerOut="noWicket"
    }

    if(plot == TRUE){
        plot.title = paste("Wicket kind taken by bowlers -",main," Vs ",opposition,"(all matches)",sep="")
        ggplot(data=r,aes(x=wicketKind,y=runs,fill=factor(wicketKind))) +
            facet_wrap( ~ bowler,scales = "fixed", ncol=8) +
            geom_bar(stat="identity") +
            xlab("Wicket kind") + ylab("Runs conceded") +
            ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    else{
        r
    }


}
