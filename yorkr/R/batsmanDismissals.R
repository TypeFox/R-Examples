##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: batsmanDismissals
# This function computes and plots the dismissal types for the batsman
#
###########################################################################################
#' @title
#' Dismissal type of batsmen
#'
#' @description
#' This function computes and plots the type of dismissals of the
#' the batsman
#' @usage
#' batsmanDismissals(df,name="A Leg Glance")
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
#' batsmanDismissals(kohli,"Kohli")
#' }
#' @seealso
#' \code{\link{batsmanFoursSixes}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#'
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#' @import yaml
#' @import rpart.plot
#' @importFrom gridExtra grid.arrange
#' @importFrom stats complete.cases loess
#' @importFrom utils head
#' @importFrom graphics legend lines plot
#'
#'
#' @export
#'

batsmanDismissals <- function(df,name="A Leg Glance"){
    batsman <- wicketKind <-dismissal <- NULL
    DismissalType <- NULL
    b <-select(df,batsman,wicketKind)

    c <- summarise(group_by(b,batsman,wicketKind),dismissal=n())
    d <- mutate(c,wicketKind=paste(wicketKind,"(",dismissal,")",sep=""))
    names(d) <-c("batsman","DismissalType","dismissal")
    plot.title = paste(name,"- Dismissals")
    ggplot(d, aes(x="", y=dismissal, fill=DismissalType))+
        geom_bar(width=1,stat = "identity")+
        coord_polar("y", start=0) +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))


}
