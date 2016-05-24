##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 20 Mar 2016
# Function: getMatchDetails
# This function loads the data for a specified match given the teams in the match
# and the date the match was played
#
###########################################################################################
#' @title
#' Get  match details of 2 countries
#'
#' @description
#' This function  gets the details of a matc palyed between 2 countries from the saved RData files
#' and returns a dataframe
#'
#' @usage
#' getMatchDetails(team1,team2,date,dir=".")
#'
#' @param team1
#' The 1st team in the match
#'
#' @param team2
#' The 2nd team in the matcj
#'
#' @param date
#' The date on which the match was played
#'
#' @param dir
#' The source directory of the RData files with all matches
#'
#' @return match
#' The dataframe of the match
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
#' # convertAllYaml2RDataframes() & convertYaml2RDataframe convert yaml files
#' # to data frame and store as RData
#' # We have to point to this directory for the call below
#' a <- getMatchDetails("England","Pakistan","2006-09-05",dir="../data")
#'
#' # Use this to create a classification tree of deliveries to wickets
#' bowlerWktsPredict(jadeja1,"RA Jadeja")
#' }
#'
#' @seealso
#' \code{\link{getBatsmanDetails}}\cr
#' \code{\link{getBowlerWicketDetails}}\cr
#' \code{\link{getTeamBattingDetails}}\cr
#' \code{\link{getTeamBowlingDetails}}\cr
#'
#' @export
#'
getMatchDetails <- function(team1,team2,date,dir="."){
    overs <- NULL
    match <- NULL
    # Create 2 filenames with both combinations of team1 and team2
    d1 <- paste(team1,"-",team2,"-",date,".RData",sep="")
    d2 <- paste(team2,"-",team1,"-",date,".RData",sep="")
    path1=paste(dir,"/",d1,sep="")
    path2=paste(dir,"/",d2,sep="")
    if(file.exists(path1)){
        load(path1)
        match <- overs
    } else if(file.exists(path2)){
        load(path2)
        match <- overs
    }else {
        cat("Match file not found at",dir, "\n")
    }

}
