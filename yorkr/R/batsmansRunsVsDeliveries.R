##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: batsmanRunsVsDeliveries
# This function plots the runs scored vs deliveries faced
#
###########################################################################################
#' @title
#' Runs versus deliveries faced
#'
#' @description
#' This function plots the runs scored and the deliveries required. A  regression
#' smoothing function is used to fit the points
#'
#' @usage
#' batsmanRunsVsDeliveries(df, name= "A Late Cut")
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
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' #Get the data frame for Kohli
#' kohli <- getBatsmanDetails(team="India",name="Kohli",dir=pathToFile)
#' batsmanRunsVsDeliveries(kohli,"Kohli")
#' }
#'
#' @seealso
#' \code{\link{batsmanFoursSixes}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#'
#' @export
#'

batsmanRunsVsDeliveries <- function(df,name= "A Late Cut"){
    batsman = runs  = ballsPlayed= NULL

    plot.title = paste(name,"- Runs vs balls faced")
    ggplot(df,aes(x=ballsPlayed,y=runs)) +
        geom_point(size=2) + geom_smooth() +
        xlab("Deliveries faced") + ylab("Runs") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))

}
