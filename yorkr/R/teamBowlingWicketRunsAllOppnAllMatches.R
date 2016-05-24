##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBowlingWicketRunsAllOppnAllMatches
# This function computes the wickets taken and the runs conceded  against all  opposition
#
###########################################################################################
#' @title
#' Team bowling wicket runs all matches against all oppositions
#'
#' @description
#' This function computes the number of wickets and runs conceded by bowlers in all matches against
#' all oppositions. The user can chose to plot or return a data frame
#'
#' @usage
#' teamBowlingWicketRunsAllOppnAllMatches(matches,t1,t2="All",plot=TRUE)
#'
#' @param matches
#' The matches of the team against all oppositions and all matches
#'
#' @param t1
#' Team for which bowling performance is required
#'
#' @param t2
#' t2=All gives the performance of the team against all opponents. Giving a opposing team (Australia, India
#' ) will give the performance against this  team
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
#' #Get all matches between India  and other opposition
#' matches <-getAllMatchesAllOpposition("India",dir="../data/",save=TRUE)
#'
#' # Or load directly from saved file
#' # load("allMatchesAllOpposition-India.RData")
#'
#' teamBowlingWicketRunsAllOppnAllMatches(matches,t1="India",t2="All",plot=TRUE)
#' m <-teamBowlingWicketRunsAllOppnAllMatches(matches,t1="India",t2="All",plot=FALSE)
#' }
#'
#' @seealso
#' \code{\link{teamBowlingScorecardAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesPlot}}\cr
#'
#' @export
#'
teamBowlingWicketRunsAllOppnAllMatches <- function(matches,t1,t2="All",plot=TRUE){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=NULL
    over=wickets=NULL

    a <- NULL
    if(t2 == "All"){
        a <-filter(matches,team==t1)
    } else {
        a <-filter(matches,team==t2)
    }

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

    #Compute number of wickets (remove nobody)
    c <- b %>%
        select(bowler,wicketKind,wicketPlayerOut) %>%
        filter(wicketPlayerOut != "nobody")

    # Count wickets by bowlers
    d <- summarise(group_by(c,bowler),wickets=length(wicketPlayerOut))

    # Calculate runs
    e <- summarise(group_by(b,bowler,over),sum(runs,wides,noballs))
    names(e) <- c("bowler","over","runs")


    #Compute total runs conceded (runs_wides+noballs)
    f <- summarize(group_by(e,bowler),runsConceded=sum(runs))

    # Join the runs conceded with the wickets taken
    g <- full_join(f,d,by="bowler")

    # Set the NAs (0 wickets) to 0
    if(sum(is.na(g$wickets)) != 0){
        g[is.na(g$wickets),]$wickets=0
    }

    # Pick the top 10 bowlers
    h <- arrange(g,desc(wickets))
    k <- h[1:10,]

    if(plot==TRUE){
        plot.title <- paste(t1,"vs",t2,"wicket Runs of bowlers")
        ggplot(data=k,aes(x=factor(wickets),y=runsConceded,fill=factor(wickets))) +
            facet_grid( ~ bowler) + geom_bar(stat="identity") +
            xlab("Number of wickets") + ylab('Runs conceded') +
            ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }else{
        k
    }

}
