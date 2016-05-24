##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBowlersVsBatsmenAllOppnAllMatchesMain
# This function computes the performance of bowlers of team against all opposition in all matches
# This function returns a dataframe
#
###########################################################################################
#' @title
#' Compute team bowlers vs batsmen all opposition  all matches
#'
#' @description
#' This function computes performance of bowlers of a team against all opposition in all matches
#'
#' @usage
#' teamBowlersVsBatsmenAllOppnAllMatchesMain(matches,theTeam,rank=0)
#'
#' @param matches
#' the data frame of all matches between a team and aall opposition  and all obtained with
#' the call getAllMatchesAllOpposition()
#'
#' @param theTeam
#' The team against which the performance is requires
#'
#' @param rank
#' When the rank is 0 then the performance of all the bowlers is displayed. If rank=n (1,2,3 ..) then
#' the performance of that bowler is given
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
#' # Get all matches between India and all oppostions
#' matches <-getAllMatchesAllOpposition("India",dir="../data/",save=TRUE)
#'
#' # You could also load directly from the saved file
#' #load("allMatchesAllOpposition-India.RData")
#' # The call below gives the best bowlers of India
#' teamBowlersVsBatsmenAllOppnAllMatchesMain(matches,theTeam="India",rank=0)
#'
#' # The call with rank=1 gives the performance of the 'India' bowler with rank=1
#'  teamBowlersVsBatsmenAllOppnAllMatchesMain(matches,theTeam="India",rank=1)
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
teamBowlersVsBatsmenAllOppnAllMatchesMain <- function(matches,theTeam,rank=0) {
    team=bowler=batsman=NULL
    runs=over=runsConceded=NULL

    a <-filter(matches,team !=theTeam)

    b <-summarise(group_by(a,bowler,batsman),sum(runs))
    names(b) <- c("bowler","batsman","runsConceded")
    # Compute total runs conceded
    c <- summarise(group_by(b,bowler),runs=sum(runsConceded))
    # Sort by descneding
    d <- arrange(c,desc(runs))


    # Initialise to NULL
    f <- NULL
    if(rank == 0){
        f <- head(d,10)
    } else { # display dispRows for selected bowler with rank
       # Pick the chosen bowler
        bwlr <- d[rank,]

        f <- filter(b,bowler==bwlr$bowler)
        f <- arrange(f,desc(runsConceded))
    }
    f


}
