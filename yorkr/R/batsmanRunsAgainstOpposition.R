##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 26 Mar 2016
# Function: batsmanRunsAgainstOpposition
# This function computes and plots the runs scored by the batsman against different oppositions
#
###########################################################################################
#' @title
#' Batsman runs against different oppositions
#'
#' @description
#' This function computes and plots the mean runs scored by the batsman against different
#' oppositions
#' @usage
#' batsmanRunsAgainstOpposition(df, name= "A Leg Glance")
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
#' batsmanRunsAgainstOpposition(kohli,"Kohli")
#' }
#'
#' @seealso
#' \code{\link{batsmanFoursSixes}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#' \code{\link{batsmanRunsPredict}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#'
#' @export
#'

batsmanRunsAgainstOpposition <- function(df,name= "A Leg Glance"){
    batsman = runs = opposition = meanRuns =  NULL
    b <- select(df,batsman,runs,opposition)
    c <-b[complete.cases(b),]
    d <- summarise(group_by(c,opposition),meanRuns=mean(runs),numMatches=n())
    plot.title = paste(name,"- Runs against opposition")
    ggplot(d, aes(x=opposition, y=meanRuns, fill=opposition))+
        geom_bar(stat = "identity",position="dodge") +
        xlab("Opposition") + ylab("Runs") +
        geom_hline(aes(yintercept=50))+
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
