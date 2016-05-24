##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 22 Mar 2016
# Function: teamBatsmenVsBowlersOppnAllMatches
# This function computes the best performing batsman against an opposition's bowlers
# in all matches with this team.The top 5 batsman are displayed by default
#
#
###########################################################################################
#' @title
#' Team batsmen vs bowlers all matches of an opposition
#'
#' @description
#' This function computes the performance of batsmen against the  bowlers of an oppositions in all matches
#'
#' @usage
#' teamBatsmenVsBowlersOppnAllMatches(matches,main,opposition,plot=TRUE,top=5)
#'
#' @param matches
#' All the matches of the team against one specific opposition
#'
#' @param main
#' The team for which the the batting partnerships are sought
#'
#' @param opposition
#' The opposition team
#'
#' @param plot
#' If plot=TRUE then a plot will be displayed else a data frame will be returned
#'
#' @param top
#' The number of players to be plotted or returned as a dataframe. The default is 5
#'
#'
#' @return None or dataframe
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
#' # Get all matches for team India against an opposition
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#'
#' # Get the performance of India batsman against Australia in all matches
#' teamBatsmenVsBowlersOppnAllMatches(a,"India","Australia")
#'
#' # Display top 3
#' teamBatsmanVsBowlersOppnAllMatches(a,"Australia","India",top=3)
#'
#' # Get top 10 and do not plot
#' n <- teamBatsmenVsBowlersOppnAllMatches(a,"Australia","India",top=10,plot=FALSE)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenVsBowlersOppnAllMatches}}\cr
#'
#' @export
#'
teamBatsmenVsBowlersOppnAllMatches <- function(matches,main,opposition,plot=TRUE,top=5){
    team=batsman=bowler=runs=runsScored=NULL
    a <-filter(matches,team==main)
    b <-summarise(group_by(a,batsman,bowler),sum(runs))
    names(b) <- c("batsman","bowler","runs")
    c <- summarise(b,runsScored=sum(runs))
    d <- arrange(c,desc(runsScored))

    # Pick 9 highest run givers
    d <- head(d,top)

    batsmen <- as.character(d$batsman)
    e <- NULL
    for(i in 1:length(batsmen)){
        f <- filter(b,batsman==batsmen[i])
        e <- rbind(e,f)

    }
    if(plot == TRUE){
        plot.title = paste("Batsmen vs bowlers -",main," Vs ",opposition,"(all matches)",sep="")
        ggplot(data=e,aes(x=bowler,y=runs,fill=factor(bowler))) +
            facet_grid(~ batsman) + geom_bar(stat="identity") +
            xlab("Bowler") + ylab("Runs Scored") +
            ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    else{
        e
    }
}

