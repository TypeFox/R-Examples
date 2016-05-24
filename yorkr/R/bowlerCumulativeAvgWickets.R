##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 14 Apr 2016
# Function: bowlerCumulativeAvgWickets
# This function computes and plots the cumulative average wickets of a bowler
#
###########################################################################################
#' @title
#' Bowler's cumulative average wickets
#'
#' @description
#' This function computes and plots the cumulative average wickets of a bowler
#'
#' @usage
#' bowlerCumulativeAvgWickets(df,name)
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
#' bowlerCumulativeAvgWickets(jadeja,"RA Jadeja")
#' }
#' @seealso
#' \code{\link{batsmanCumulativeAverageRuns}}
#' \code{\link{bowlerCumulativeAvgEconRate}}
#' \code{\link{batsmanCumulativeStrikeRate}}
#' \code{\link{batsmanRunsVsStrikeRate}}
#' \code{\link{batsmanRunsPredict}}
#'
#' @export
#'
bowlerCumulativeAvgWickets <- function(df,name){
    wickets=cs=no=NULL
    b <- select(df,wickets)
    b$no<-seq.int(nrow(b))
    c <- select(b,no,wickets)

    d <- mutate(c,cs=cumsum(wickets)/no)
    plot.title= paste(name,"- Cumulative avg wkts vs No innings")
    ggplot(d) + geom_line(aes(x=no,y=cs),col="blue") +
        xlab("No of innings") + ylab("Cumulative Avg. wickets") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))
}
