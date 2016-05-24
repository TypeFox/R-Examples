##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 22 Mar 2016
# Function: teamBatsmenPartnershipOppnAllMatchesChart
# This function computes the batting partnership of a team in all matches against
# an opposition. The report generated can be detailed or a summary
#
###########################################################################################
#' @title
#' Plot of team partnership all matches against an opposition
#'
#' @description
#' This function plots  the batting partnership of a team againt all oppositions in all matches
#' This function also returns a  dataframe  with the batting partnerships
#'
#' @usage
#' teamBatsmenPartnershipOppnAllMatchesChart(matches,main,opposition, plot=TRUE)
#'
#' @param matches
#' All the matches of the team against all oppositions
#'
#' @param main
#' The main team for which the the batting partnerships are sought
#'
#' @param opposition
#' The opposition team for which the the batting partnerships are sought
#'
#' @param plot
#' Whether the partnerships have top be rendered as a plot. If plot=FALSE the data frame is returned
#'
#' @return None or partnerships
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
#' # Get all matches for team India against all oppositions
#' d <- teamBatsmenVsBowlersAllOppnAllMatchesRept(matches,"India",rank=1,dispRows=50)
#'  #Plot the partnerships
#' teamBatsmenVsBowlersAllOppnAllMatchesPlot(d)
#'
#' #Do not plot but get the dataframe
#' e <- teamBatsmenVsBowlersAllOppnAllMatchesPlot(d,plot=FALSE)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenVsBowlersOppnAllMatches}}\cr
#'
#' @export
#'
teamBatsmenPartnershipOppnAllMatchesChart <- function(matches,main,opposition,plot=TRUE){
    team=batsman=nonStriker=runs=partnershipRuns=totalRuns=NULL
    a <-filter(matches,team==main)
    #Get partnerships
    df <- data.frame(summarise(group_by(a,batsman,nonStriker),sum(runs)))
    names(df) <- c("batsman","nonStriker","runs")

    df <- arrange(df,desc(runs))

    plot.title = paste("Batting partnership -",main," Vs ",opposition,"(all matches)",sep="")
    if(plot==TRUE){
    # Plot the data
    ggplot(data=df,aes(x=batsman,y=runs,fill=nonStriker))+
        geom_bar(data=df,stat="identity") +
        xlab("Batsman") + ylab("Partnership runs") +
        ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } else{
        df
    }


}
