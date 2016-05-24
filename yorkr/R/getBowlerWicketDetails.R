##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: getBowlerWicketDetails
# This function gets the bowling details of a bowler
#
###########################################################################################
#' @title
#' Get the bowling details of a bowler
#'
#' @description
#' This function gets the bowling of a bowler (overs,maidens,runs,wickets,venue, opposition)
#'
#' @usage
#' getBowlerWicketDetails(team,name,dir=".")
#'
#' @param team
#' The team to which the bowler belongs
#'
#' @param name
#' The name of the bowler
#'
#' @param dir
#' The source directory of the data
#'
#' @return dataframe
#' The dataframe of bowling performance
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
#'
#' @examples
#' \dontrun{
#' # Get the bowling details of bowlers of a team e.g. India. This is saved as a dataframe
#' c <- getTeamBowlingDetails("India",dir="../data",save=TRUE)
#' #Get the bowler details from this overall data frame
#'
#' jadeja <- getBowlerWicketDetails(team="India",name="Jadeja",dir=".")
#'
#' # The dataframe from the above call is used in many functions
#' #bowlerMeanEconomyRate(jadeja,"RA Jadeja")
#' }
#'
#' @seealso
#' \code{\link{bowlerMovingAverage}}\cr
#' \code{\link{getTeamBowlingDetails}}\cr
#' \code{\link{bowlerMeanRunsConceded}}\cr
#' \code{\link{teamBowlersWicketRunsOppnAllMatches}}\cr
#'
#' @export
#'
getBowlerWicketDetails <- function(team,name,dir="."){
    bowlingDetails=bowler=wicketPlayerOut=overs=maidens=NULL
    runs=economyRate=opposition=wickets=venue=NULL
    fl <- paste(dir,"/",team,"-BowlingDetails.RData",sep="")
    load(fl)
    details <- bowlingDetails
    bowlerDetails <- filter(details,grepl(name,bowler))
    bowlerDetails <- arrange(bowlerDetails,date)

    # Count wickets taken
    # Replace nobody by 0 wickets
    a <- filter(bowlerDetails,wicketPlayerOut == "nobody")
    a$wickets <- 0
    a1 <- c <- select(a,bowler,overs,maidens,runs,economyRate,date,
                      opposition,wickets,venue)

    # Get rows which have wickets
    b <- filter(bowlerDetails,wicketPlayerOut != "nobody")
    c <- select(b,bowler,overs,maidens,runs,economyRate,date,opposition,venue)
    # Count wickets
    d <- summarise(group_by(b,date),wickets=length(unique(wicketPlayerOut)))
    # Join tables
    e <- full_join(c,d,by="date")

    f <- rbind(a1,e)
    f <- select(f,bowler,overs,maidens,runs,wickets,economyRate,date,opposition,venue)
    g <- arrange(f,date)
    g


}
