##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: batsmanMovingAverage
# This function computes and plots the moving average  the batsman
#
###########################################################################################
#' @title
#' Moving average of batsman
#'
#' @description
#' This function plots the runs scored by the batsman over the career as a time
#' series. A loess regression line is plotted on the moving average of the batsman
#' the batsman
#'
#' @usage
#' batsmanMovingAverage(df, name= "A Leg Glance")
#'
#' @param df
#' Data frame
#'
#' @param name
#' Name of batsman
#'
#' @return None
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
#' #Get the data frame for Kohli
#' kohli <- getBatsmanDetails(team="India",name="Kohli",dir=pathToFile)
#' batsmanMovingAverage(kohli,"Kohli")
#' }
#' @seealso
#' \code{\link{batsmanDismissals}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#' \code{\link{batsmanRunsPredict}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#'
#' @export
#'
batsmanMovingAverage <- function(df,name = "A Leg Glance"){
    batsman = runs = NULL
    b <- select(df,batsman,runs,date)

    plot.title = paste(name,"- Moving average of runs in career")
    ggplot(b) + geom_line(aes(x=date, y=runs),colour="darkgrey") +
        geom_smooth(aes(x=date, y=runs)) +
        xlab("Date") + ylab("Runs") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
