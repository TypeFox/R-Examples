##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 16 Apr 2016
# Function: batsmanRunsVenue
# This function computes and plots the runs scored by the batsman at different venues
#
###########################################################################################
#' @title
#' Batsman runs at different venues
#'
#' @description
#' This function computes and plots the mean runs scored by the batsman at different
#' venues of the world
#' @usage
#' batsmanRunsVenue(df, name= "A Leg Glance")
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
#' batsmanRunsVenue(kohli,"Kohli")
#' }
#'
#' @seealso
#' \code{\link{batsmanFoursSixes}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#' \code{\link{batsmanRunsPredict}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{batsmanRunsAgainstOpposition}}\cr
#'
#' @export
#'

batsmanRunsVenue <- function(df,name= "A Leg Glance"){
    batsman = runs = venue = numMatches = meanRuns = NULL
    b <- select(df,batsman,runs,venue)
    c <- summarise(group_by(b,venue),meanRuns=mean(runs),numMatches=n())
    d <- mutate(c,venue=paste(venue,"(",numMatches,")",sep=""))
    e <- arrange(d,desc(numMatches))

    # Select only available rows or 25 rows
    sz <- dim(e)
    if(sz[1] > 25){
        f <- e[1:25,]
    } else{
        f <- e[1:sz[1],]
    }
    plot.title = paste(name,"- Mean runs at venue")
    ggplot(f, aes(x=venue, y=meanRuns, fill=venue))+
        geom_bar(stat = "identity",position="dodge") +
        geom_hline(aes(yintercept=50))+
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

}
