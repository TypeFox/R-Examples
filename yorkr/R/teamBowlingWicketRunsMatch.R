##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 21 Mar 2016
# Function: teamBowlingWicketRunsMatch
# This function computes the performance of bowlers of team with the runs conceded
#
#
###########################################################################################
#' @title
#' Team bowling wickets runs conceded in match
#'
#' @description
#' This function computes returns the wickets taken and runs conceded bowlers in a match between 2 teams.
#' The user can choose to plot or return a dataframe
#'
#' @usage
#' teamBowlingWicketRunsMatch(match,theTeam,plot=TRUE)
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
#' teamBowlingWicketRunsMatch(a,"England",plot=FALSE)
#' teamBowlingWicketRunsMatch(a,"Pakistan")
#' }
#'
#' @seealso
#' \code{\link{teamBowlingWicketMatch}}\cr
#' \code{\link{teamBowlingWicketRunsMatch}}\cr
#' \code{\link{teamBowlersVsBatsmenMatch}}\cr
#'
#' @export
#'
teamBowlingWicketRunsMatch <- function(match,theTeam,plot=TRUE){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=over=wickets=NULL
    # The performance of bowlers of the team is got when the other side is batting. Hence '!-"
    # Filter the bowler's performance
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

    l$wickets = as.character(l$wickets)
    # Set NAs to 0
    if(sum(is.na(l$wickets)) != 0){
        l[is.na(l$wickets),]$wickets=0
    }

    # Plot or ourput data frame
    if(plot == TRUE){
        plot.title <- paste(theTeam,"Number of wickets vs Runs conceded by bowlers")
        ggplot(data=l,aes(x=factor(wickets),y=runs,fill=factor(wickets))) +
            facet_grid(. ~ bowler,scales = "free_x", space = "free_x") +
            geom_bar(stat="identity") +
            xlab("Number of wickets") + ylab("Total runs conceded") +
            ggtitle(bquote(atop(.(plot.title),
                                    atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    else {
        l
    }

}
