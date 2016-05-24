##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 14 Apr 2016
# Function: batsmanCumulativeStrikeRate
# This function computes and plots the cumulative average strike rate of a batsman
#
###########################################################################################
#' @title
#' Batsman's cumulative average strike rate
#'
#' @description
#' This function computes and plots the cumulative average strike rate  of a batsman
#'
#' @usage
#' batsmanCumulativeStrikeRate(df,name= "A Leg Glance")
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
#' batsmanCumulativeStrikeRate(kohli,"Kohli")
#' }
#' @seealso
#' \code{\link{batsmanCumulativeAverageRuns}}
#' \code{\link{bowlerCumulativeAvgEconRate}}
#' \code{\link{bowlerCumulativeAvgWickets}}
#' \code{\link{batsmanRunsVsStrikeRate}}
#' \code{\link{batsmanRunsPredict}}
#'
#' @export
#'
batsmanCumulativeStrikeRate <- function(df,name= "A Leg Glance"){
    strikeRate=cs=no=NULL
    b <- select(df,strikeRate)
    b$no<-seq.int(nrow(b))
    c <- select(b,no,strikeRate)

    d <- mutate(c,cs=cumsum(strikeRate)/no)
    plot.title= paste(name,"- Cumulative Strike Rate vs No of innings")
    ggplot(d) + geom_line(aes(x=no,y=cs),col="blue") +
        xlab("No of innings") + ylab("Cumulative Avg. Strike Rates") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
