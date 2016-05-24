##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 21 Mar 2016
# Function: teamBowlingWicketMatch
# This function computes the performance of bowlers and the wickets taken
# The user has a choice of either taking the output as a plot or as a dataframe
#
###########################################################################################
#' @title
#' Compute and plot wickets  by bowlers in  match
#'
#' @description
#' This function computes returns the wickets taken bowlers in a match between 2 teams
#'
#' @usage
#' teamBowlingWicketMatch(match,theTeam,plot=TRUE)
#'
#' @param match
#' The match between the teams
#'
#' @param theTeam
#' Team for which bowling performance is required
#'
#' @param plot
#' If plot= TRUE the dataframe will be plotted else a data frame will be returned
#'
#' @return None or data fame
#' A data frame with the bowling performance in alll matches against all oppositions
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
#' #Get the match details
#' a <- getMatchDetails("England","Pakistan","2006-09-05",dir="../temp")
#'
#' teamBowlingWicketMatch(a,"England",plot=FALSE)
#' teamBowlingWicketMatch(a,"Pakistan")
#' }
#'
#' @seealso
#' \code{\link{teamBowlingWicketMatch}}\cr
#' \code{\link{teamBowlingWicketRunsMatch}}\cr
#' \code{\link{teamBowlersVsBatsmenMatch}}\cr
#'
#' @export
#'
teamBowlingWicketMatch <- function(match,theTeam,plot=TRUE){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=over=NULL
    # The bowlers performance of the team is got when the other side is batting. Hence '!-"
    # Filter the data frame
    a <-filter(match,team!=theTeam)
    # Compute the maidens,runs conceded and overs for the bowlers
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


    #Compute number of wickets
    h <- b %>%
        select(bowler,wicketKind,wicketPlayerOut) %>%
        filter(wicketPlayerOut != "nobody")
    i <- summarise(group_by(h,bowler),wickets=length(unique(wicketPlayerOut)))

    j <- full_join(h,e,by="bowler")

    # Make as.character to assign values
    j$wicketKind = as.character(j$wicketKind)
    j$wicketPlayerOut = as.character(j$wicketPlayerOut)
    # Set NAs to 0
    if(sum(is.na(j$wicketKind) !=0)){
        j[is.na(j$wicketKind),]$wicketKind="noWicket"
    }
    if(sum(is.na(j$wicketPlayerOut) != 0)){
        j[is.na(j$wicketPlayerOut),]$wicketPlayerOut="noWicket"
    }

    if(plot == TRUE){
        plot.title <- paste(theTeam,"Wickets taken vs Runs conceded by bowlers")
        ggplot(data=j,aes(x=wicketPlayerOut,y=runs,fill=factor(wicketPlayerOut))) +
            facet_grid(. ~ bowler,scales = "free_x", space = "free_x") +
            geom_bar(stat="identity") +
            xlab("Batsman out") + ylab("Total runs conceded") +
            ggtitle(bquote(atop(.(plot.title),
                                    atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    else{
        j
    }


}
