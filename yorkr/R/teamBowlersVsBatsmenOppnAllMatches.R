##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 23 Mar 2016
# Function: teamBowlersVsBatsmenOppnAllMatches
# This function computes the performance of the bowlers and the runs conceded and the batsman
# who scored most
#
#
###########################################################################################
#' @title
#' Team bowlers vs batsmen against an opposition in all matches
#'
#' @description
#' This function computes performance of bowlers of a team against an opposition in all matches
#' against the opposition
#'
#' @usage
#' teamBowlersVsBatsmenOppnAllMatches(matches,main,opposition,plot=TRUE,top=5)
#'
#' @param matches
#' The data frame of all matches between a team the opposition. This dataframe can be obtained with
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#'
#' @param main
#' The main team against which the performance is requires
#'
#' @param opposition
#' The opposition team against which the performance is require
#'
#' @param plot
#' If true plot else return dataframe
#'
#' @param top
#' The number of rows to be returned. 5 by default
#'
#' @return dataframe
#' The dataframe with all performances
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
#' # Get all matches between India and Australia
#' matches <- getAllMatchesBetweenTeams("Australia","India",dir="../data")
#'
#' #  Plot the performance of top 5 Indian bowlers against Australia
#' teamBowlersVsBatsmanOppnAllMatches(matches,'India',"Australia",top=5)
#'
#' # Plot the performance of top 3 Australian bowlers against India
#' teamBowlersVsBatsmenOppnAllMatches(matches,"Australia","India",top=3)
#'
#' # Get the top 5 bowlers of Australia. Do not plot but get as a dataframe
#' teamBowlersVsBatsmenOppnAllMatches(matches,"Australia","India",plot=FALSE)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesRept}}\cr
#' \code{\link{teamBowlersVsBatsmenAllOppnAllMatchesPlot}}\cr
#'
#' @export
#'
teamBowlersVsBatsmenOppnAllMatches <- function(matches,main,opposition,plot=TRUE,top=5){
    noBalls=wides=team=runs=bowler=wicketKind=wicketPlayerOut=NULL
    team=bowler=ball=wides=noballs=runsConceded=overs=batsman=NULL
    a <-filter(matches,team != main)

    b <-summarise(group_by(a,bowler,batsman),sum(runs))
    names(b) <- c("bowler","batsman","runsConceded")
    # Compute total runs conceded
    c <- summarise(group_by(b,bowler),runs=sum(runsConceded))
    # Sort by descneding
    d <- arrange(c,desc(runs))

    # Pick 5 highest run givers
    d <- head(d,top)

    bowlers <- as.character(d$bowler)
    e <- NULL
    for(i in 1:length(bowlers)){
        f <- filter(b,bowler==bowlers[i])
        e <- rbind(e,f)

    }
    names(e) <- c("bowler","batsman","runsConceded")

    if( plot == TRUE){
        plot.title = paste("Bowlers vs batsmen -",main," Vs ",opposition,"(all matches)",sep="")
        ggplot(data=e,aes(x=batsman,y=runsConceded,fill=factor(batsman))) +
            facet_grid(. ~ bowler) + geom_bar(stat="identity") +
            #facet_wrap( ~ bowler,scales = "free", ncol=3,drop=TRUE) + #Does not work.Check!
            xlab("Batsman") + ylab("Runs conceded") +
            ggtitle(bquote(atop(.(plot.title),
                                    atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } else{
        e
    }
}
