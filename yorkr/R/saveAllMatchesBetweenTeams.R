##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 30 Mar 2016
# Function: saveAllMatchesBetweenTeams
# This function saves all matches between 2 teams as a single dataframe
##################################################################################
#' @title
#' Saves all matches between 2 teams as  dataframe
#'
#' @description
#' This function saves all matches between 2 teams as a single dataframe in the
#' current directory
#'
#' @usage
#' saveAllMatchesBetweenTeams()
#'
#' @return None
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#'
#'
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' saveAllMatchesBetweenTeams
#' }
#' @seealso
#' \code{\link{batsmanDismissals}}\cr
#' \code{\link{batsmanRunsVsDeliveries}}\cr
#' \code{\link{batsmanRunsVsStrikeRate}}\cr
#' \code{\link{getAllMatchesAllOpposition}}\cr
#' \code{\link{getAllMatchesBetweenTeams}}\cr
#'
#' @export
#'

saveAllMatchesBetweenTeams <- function(){

    teams <-c("Australia","India","Pakistan","West Indies", 'Sri Lanka',
              "England", "Bangladesh","Netherlands","Scotland", "Afghanistan",
              "Zimbabwe","Ireland","New Zealand","South Africa","Canada",
              "Bermuda","Kenya")

    matches <- NULL
    #Create all combinations of teams
    for(i in seq_along(teams)){
        for(j in seq_along(teams)){
            if(teams[i] != teams[j]){
                cat("Team1=",teams[i],"Team2=",teams[j],"\n")
                tryCatch(matches <- getAllMatchesBetweenTeams(teams[i],teams[j],dir="../data",save=TRUE),
                         error = function(e) {
                             print("No matches")

                         }
                )
            }
        }
        matches <- NULL
    }
}

