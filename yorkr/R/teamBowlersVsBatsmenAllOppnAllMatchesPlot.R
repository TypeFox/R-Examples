##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBowlersVsBatsmenAllOppnAllMatchesPlot
# This function computes the performance of bowlers against batsman of opposition
#
###########################################################################################
#' @title
#' Plot bowlers vs batsmen against all opposition all matches
#'
#' @description
#' This function computes performance of bowlers of a team against all opposition in all matches
#'
#' @usage
#' teamBowlersVsBatsmenAllOppnAllMatchesPlot(bowlerDF,t1,t2)
#'
#' @param bowlerDF
#' The data frame of the bowler whose performance is required
#'
#' @param t1
#' The team against to which the player belong
#'
#' @param t2
#' The opposing team
#'
#'
#' @return none
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
#' # Get all matches between India and all oppostions
#' matches <-getAllMatchesAllOpposition("India",dir="../data/",save=TRUE)
#'
#' #Get the details of the bowler with the specified rank as a dataframe
#' df <- teamBowlersVsBatsmenAllOppnAllMatchesRept(matches,theTeam="India",rank=1)
#' #Plot this
#' teamBowlersVsBatsmenAllOppnAllMatchesPlot(df,"India","India")
#'
#' df <- teamBowlersVsBatsmenAllOppnAllMatchesRept(matches,theTeam="England",rank=1)
#' teamBowlersVsBatsmenAllOppnAllMatchesPlot(df,"India","England")
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesMain}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesRept}}\cr
#' @export
#'
teamBowlersVsBatsmenAllOppnAllMatchesPlot <- function(bowlerDF,t1,t2){
    batsman=runsConceded=team=NULL
    bwlr <- bowlerDF$bowler
    if(t2 != "India"){
        plot.title <- paste(bwlr,"-Performance against",t2,"batsmen")
        print("aa")
    }else{
        plot.title <- paste(bwlr,"-Performance against all batsmen")
    }
    ggplot(data=bowlerDF,aes(x=batsman,y=runsConceded,fill=factor(batsman))) +
        facet_grid(. ~ bowler) + geom_bar(stat="identity") +
        xlab("Batsman") + ylab("Runs conceded") +
        ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
