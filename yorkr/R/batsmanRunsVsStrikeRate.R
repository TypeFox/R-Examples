##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: batsmanRunsVsStrikeRate
# This function plots the runs scored versus the strike rate
###########################################################################################
#' @title
#' Batsman runs versus strike rate
#'
#' @description
#' This function plots the runs scored by the batsman  and the runs scored
#' by the batsman. A loess line is fitted over the points
#'
#' @usage
#' batsmanRunsVsStrikeRate(df, name= "A Late Cut")
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
#' batsmanRunsVsStrikeRate(kohli,"Kohli")
#' }
#'
#' @seealso
#' \code{\link{batsmanDismissals}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#' \code{\link{batsmanRunsPredict}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#'
#' @export
#'
batsmanRunsVsStrikeRate <- function(df,name= "A Late Cut"){
    batsman = runs = strikeRate = NULL
    b <- select(df,batsman,runs,strikeRate)

    plot.title = paste(name,"- Runs vs Strike Rate")
    ggplot(b) + geom_point(aes(x=runs, y=strikeRate),colour="darkgrey") +
        geom_smooth(aes(x=runs, y=strikeRate)) +
        xlab("Strike rate(%)") + ylab("Runs") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
