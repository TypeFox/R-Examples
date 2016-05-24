##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 21 Mar 2016
# Function: teamBatsmenPartnershipMatch
# This function computes and displays the partnership details in a match. The output
# can either be a plot or the data frame used in the plot
#
###########################################################################################
#' @title
#' Team batting partnerships of batsmen in a match
#'
#' @description
#' This function plots the partnerships of batsmen in a match against an opposition or it can return
#' the data frame
#'
#' @usage
#' teamBatsmenPartnershipMatch(match,theTeam,plot=TRUE)
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
#' @return df
#' The data frame of the batsmen partnetships
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
#' batsmenPartnershipMatch(a,"Pakistan")
#' batsmenPartnershipMatch(a,"England",plot=TRUE)
#' m <-batsmenPartnershipMatch(a,"Pakistan",plot=FALSE)
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
teamBatsmenPartnershipMatch <- function(match,theTeam,plot=TRUE){
    team=batsman=nonStriker=runs=runsScored=NULL
    a <-filter(match,team==theTeam)
    # Group batsman with non strikers and compute partnerships
    df <- data.frame(summarise(group_by(a,batsman,nonStriker),sum(runs)))
    names(df) <- c("batsman","nonStriker","runs")

    if(plot==TRUE){
        plot.title <- paste(theTeam,"Batting partnership in match")
        ggplot(data=df,aes(x=batsman,y=runs,fill=nonStriker))+
            geom_bar(data=df,stat="identity") +
            xlab("Batmen") + ylab("Runs Scored") +
            ggtitle(bquote(atop(.(plot.title),
                                    atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    else{
        # Output dataframe
        df
    }


}
