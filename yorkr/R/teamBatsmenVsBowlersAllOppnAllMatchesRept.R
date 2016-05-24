##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBatsmenVsBowlersAllOppnAllMatchesRept
# This function computes performance of batsmen/batsman against bowlers of the opposition.
# It provides the names of the bowlers against whom the batsman scored the most.
# We can the over all performance of the team or the individual performances of the batsman
# If rank=10 then the overall performance of the team is displayed
# For a rank'n' the performance of the batsman at that rank against bowlers is displayed
###########################################################################################
#' @title
#' Report of team batsmen vs bowlers in all matches all oppositions
#'
#' @description
#' This function computes the performance of batsmen against all bowlers of all oppositions in all matches
#'
#' @usage
#' teamBatsmenVsBowlersAllOppnAllMatchesRept(matches,theTeam,rank=0,dispRows=50)
#'
#' @param matches
#' All the matches of the team against all oppositions
#'
#' @param theTeam
#' The team for which the the batting partnerships are sought
#'
#' @param rank
#' if the rank=0  then the data frame returned gives a summary list of the batsmen with the highest
#' runs against bowlers. If rank=N (where N=1,2,3...) then the detailed breakup of the batsman and the runs
#' scored against each bowler is returned
#'
#' @param dispRows
#' The number of rows to be returned
#'
#' @return h
#' The data frame of the batsman and the runs against bowlers
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
#' # Get all matches for team India against all oppositions
#' m <-teamBattingScorecardAllOppnAllMatches(matches,theTeam="India")
#' # Get the  summary report
#' teamBatsmenVsBowlersAllOppnAllMatchesRept(matches,"India",rank=0)
#' #Get detailed report
#' teamBatsmenVsBowlersAllOppnAllMatchesRept(matches,"India",rank=1,dispRows=50)
#'
#' teamBatsmenVsBowlersAllOppnAllMatchesRept(matches,"Pakistan",rank=0)
#' teamBatsmenVsBowlersAllOppnAllMatchesRept(matches,"England",rank=1)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenVsBowlersOppnAllMatches}}\cr
#'
#' @export
#'
teamBatsmenVsBowlersAllOppnAllMatchesRept <- function(matches,theTeam,rank=0,dispRows=50)
{
    team=batsman=bowler=runs=runsScored=NULL
    a <-filter(matches,team==theTeam)
    b <-summarise(group_by(a,batsman,bowler),sum(runs))
    names(b) <- c("batsman","bowler","runs")

    c <- summarise(b,runsScored=sum(runs))
    d <- arrange(c,desc(runsScored))


    # If rank == 0 thne display top  batsman with best performance
    if(rank == 0){
        f <- d
    } else {
        # display dispRows for selected batsman with rank and runs scored against opposing bowlers
        bman <- d[rank,]
        f <- filter(b,batsman==bman$batsman)
        f <- arrange(f,desc(runs))
        # Output only dispRows
        f <- f[1:dispRows,]
    }
    g <- complete.cases(f)
    h <- f[g,]


}
