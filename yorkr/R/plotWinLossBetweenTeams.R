##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 15 Apr 2016
# Function: plotWinLossBetweenTeams
# This function computes and plots number of wins for each team
#
###########################################################################################
#' @title
#' Plot  wins for each team
#'
#' @description
#' This function computes and plots number of wins for each team in all their
#' encounters. The plot includes the number of  wins byteam1 each team and the matches
#' with no result
#'
#' @usage
#' plotWinLossBetweenTeams(team1,team2,dir=".")
#'
#' @param team1
#' The 1st team
#'
#' @param team2
#' The 2nd team
#'
#' @param dir
#' The source directory of teh RData files
#'
#' @return None
#'
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
#'
#' plotWinLossBetweenTeams(team1="India",team2="Australia",dir=pathToFile)
#' batsmanDismissals(kohli,"Kohli")
#' }
#' @seealso
#' \code{\link{batsmanFoursSixes}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#'
#'
#' @export
#'

plotWinLossBetweenTeams <- function(team1,team2,dir="."){
    overs=NULL
    venue=winner=result=date=NULL
    # Create 2 filenames with both combinations of team1 and team2
    d1 <- paste(team1,"-",team2,"*",sep="")
    d2 <- paste(team2,"-",team1,"*",sep="")
    path1=paste(dir,"/",d1,sep="")
    path2=paste(dir,"/",d2,sep="")
    # Capture both combinations
    fl1 <- Sys.glob(path1)
    fl2 <- Sys.glob(path2)
    fl3 <-c(fl1,fl2)

    # Create a data frame with all matches
    overs <- NULL
    w <- NULL
    for(i in 1:length(fl3)){
        load(fl3[i])
        o <- overs[1,]
        a <- select(o,date,venue,winner,result)
        w <- rbind(w,a)
        overs <- NULL
    }
    winLoss <- summarise(group_by(w,winner),count=n())

    x <- winLoss$winner=="NA"
    winLoss$winner <- as.character(winLoss$winner)
    if(sum(x) !=0) {
        winLoss[x,]$winner <-"NoResult"
    }

    plot.title <- paste("Number of wins in",team1," vs ",team2, " matches")
    ggplot(winLoss, aes(x=winner, y=count, fill=winner))+
        geom_bar(stat = "identity",position="dodge") +
        xlab("Winner") + ylab("Numer of Wins") +
        ggtitle(bquote(atop(.(plot.title),
                            atop(italic("Data source:http://cricsheet.org/"),""))))

}
