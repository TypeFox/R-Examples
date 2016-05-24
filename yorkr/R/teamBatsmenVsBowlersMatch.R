##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 21 Mar 2016
# Function: teamBatsmenVsBowlersMatch
# This function computes the performance of batsmen against different bowlers.
# The user has a choice of either taking the output as a plot or as a dataframe
#
###########################################################################################
#' @title
#' Team batsmen against bowlers in a match
#'
#' @description
#' This function plots the performance of  batsmen versus bowlers in a match  or it can return
#' the data frame
#'
#' @usage
#' teamBatsmenVsBowlersMatch(match,theTeam,plot=TRUE)
#'
#' @param match
#' The match between the teams
#'
#' @param theTeam
#' The team for which the the batting partnerships are sought
#'
#' @param plot
#' If plot=TRUE then a plot is created otherwise a data frame is returned
#'
#' @return b
#' The data frame of the batsmen vs bowlers performance
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
#' # Get athe match between England and Pakistan
#' a <- getMatchDetails("England","Pakistan","2006-09-05",dir="../temp")
#' batsmenVsBowlersMatch(a,'Pakistan',plot=TRUE)
#' teamBowlingScorecardMatch(a,'England')
#' teamBowlingWicketKindMatch(a,"England",plot=FALSE)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBattingScorecardMatch}}\cr
#'
#' @export
#'
teamBatsmenVsBowlersMatch <- function(match,theTeam,plot=TRUE)
{
    team=batsman=bowler=runs=runsConceded=NULL
    a <-filter(match,team==theTeam)
    # Summarise the performance of the batsmen against the bowlers vs total runs scored
    b <-summarise(group_by(a,batsman,bowler),sum(runs))
    names(b) <- c("batsman","bowler","runsConceded")

    if(plot == TRUE){
        plot.title <- paste(theTeam,"Batsmen vs Bowlers in Match")
        # Plot the performance of the batsmen as a facted grid
        ggplot(data=b,aes(x=bowler,y=runsConceded,fill=factor(bowler))) +
            facet_grid(~ batsman) + geom_bar(stat="identity") +
            xlab("Opposition bowlers") + ylab("Runs scored") +
            ggtitle('Batsmen vs Bowlers in Match') +
            ggtitle(bquote(atop(.(plot.title),
                                    atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    else{
        b
    }
}
