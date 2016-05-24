##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 26 Mar 2016
# Function: bowlerMovingAverage
# This function plots the moving average of wickets taken by bowler in career
#
###########################################################################################
#' @title
#' Bowler's moving average of wickets
#'
#' @description
#' This function computes and plots the wickets taken by the bowler over career. A loess
#' regression fit plots the moving average of wickets taken by bowler
#' @usage
#' bowlerMovingAverage(df, name)
#'
#' @param df
#' Data frame
#'
#' @param name
#' Name of bowler
#'
#' @return None
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' # Get the data frame for RA Jadeja
#' jadeja <- getBowlerWicketDetails(team="India",name="Jadeja",dir=pathToFile)
#' bowlerMeanRunsConceded(jadeja,"RA Jadeja")
#' }
#'
#' @seealso
#' \code{\link{bowlerMeanEconomyRate}}\cr
#' \code{\link{bowlerWicketPlot}}\cr
#' \code{\link{bowlerWicketsVenue}}\cr
#' \code{\link{bowlerMeanRunsConceded}}\cr
#'
#' @export
#'
bowlerMovingAverage <- function(df,name){
    bowler = wickets = NULL
    c <- select(df,bowler,wickets,date)

    plot.title = paste(name,"- Moving average of wickets in career")
    ggplot(c) + geom_line(aes(x=date, y=wickets),colour="darkgrey") +
        geom_smooth(aes(x=date, y=wickets)) +
        xlab("Date") + ylab("Wickets") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
