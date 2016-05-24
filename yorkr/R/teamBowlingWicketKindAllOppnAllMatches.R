##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBowlingWicketKindAllOppnAllMatches
# This function computes the wicket kind of bowlers against all  opposition
#
###########################################################################################
#' @title
#' team bowling wicket kind against all opposition all matches
#'
#' @description
#' This function computes returns kind of wickets (caught, bowled etc) of bowlers in all matches against
#' all oppositions. The user can chose to plot or return a data frame
#'
#' @usage
#' teamBowlingWicketKindAllOppnAllMatches(matches,t1,t2="All",plot=TRUE)
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
#' teamBowlingWicketKindAllOppnAllMatches(matches,t1="India",t2="All")
#' m <-teamBowlingWicketKindAllOppnAllMatches(matches,t1="India",t2="All",plot=FALSE)
#'
#' teamBowlingWicketKindAllOppnAllMatches(matches,t1="India",t2="Bangladesh")
#' teamBowlingWicketKindAllOppnAllMatches(matches,t1="India",t2="South Africa")
#' }
#'
#' @seealso
#' \code{\link{teamBowlingScorecardAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesPlot}}\cr
#'
#' @export
#'

teamBowlingWicketKindAllOppnAllMatches <- function(matches,t1,t2="All",plot=TRUE){

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

    # Arrange in descending order
    e <- arrange(d,desc(wickets))

    # Pick the top 8
    f <- e[1:8,]

    # Create a character vector of top 8 bowlers
    g <- as.character(f$bowler)


    # Select these top 8
    n <- NULL
    for(m in 1:8){
        mm <- filter(c,bowler==g[m])
        n <- rbind(n,mm)
    }

    # Summarise by the different wicket kinds for each bowler
    p <- summarise(group_by(n,bowler,wicketKind),m=n())

    if(plot==TRUE){
        plot.title <- paste(t1,"vs",t2,"wicket-kind of bowlers")
        # Plot
        ggplot(data=p,aes(x=wicketKind,y=m,fill=factor(wicketKind))) +
            facet_wrap( ~ bowler,scales = "fixed", ncol=8) +
            geom_bar(stat="identity") +
            xlab("Wicket kind") + ylab("Wickets") +
            ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }else{
        p
    }

}
