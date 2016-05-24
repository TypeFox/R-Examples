##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 14 Apr 2016
# Function: bowlerCumulativeAvgEconRate
# This function computes and plots the cumulative average economy rate s of a bowler
#
###########################################################################################
#' @title
#' Bowler's cumulative average economy rate
#'
#' @description
#' This function computes and plots the cumulative average economy rate  of a bowler
#'
#' @usage
#' bowlerCumulativeAvgEconRate(df,name)
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
#' \dontrun{)
#' #'Get the data frame for RA Jadeja
#' jadeja <- getBowlerWicketDetails(team="India",name="Jadeja",dir=pathToFile)
#' bowlerCumulativeAvgEconRate(jadeja,"RA Jadeja")
#' }
#' @seealso
#' \code{\link{batsmanCumulativeAverageRuns}}
#' \code{\link{bowlerCumulativeAvgWickets}}
#' \code{\link{batsmanCumulativeStrikeRate}}
#' \code{\link{batsmanRunsVsStrikeRate}}
#' \code{\link{batsmanRunsPredict}}
#'
#' @export
#'
bowlerCumulativeAvgEconRate <- function(df,name){
    economyRate=cs=no=NULL
    b <- select(df,economyRate)
    b$no<-seq.int(nrow(b))
    c <- select(b,no,economyRate)

    d <- mutate(c,cs=cumsum(economyRate)/no)
    plot.title= paste(name,"- Cum. avg Econ Rate vs No innings")
    ggplot(d) + geom_line(aes(x=no,y=cs),col="blue") +
        xlab("No of innings") + ylab("Cumulative Avg. Economy Rate") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
