##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 15 Apr 2016
# Function: getAllMatchesBetweenTeams
# This function gets all the data for matches palyed between teams and creates a large data
# frame. This data frame can be saved for suture use
# and the date the match was played
#
###########################################################################################
#' @title
#' Get data on all matches  between 2 opposing teams
#'
#' @description
#' This function gets all the data on matches between opposing teams for e.g India-Pakistan,
#' Englad-Australia, South Africa- Sri Lanka etc. It constructs a huge dataframe of all these
#' matches. This can be saved by the user which can be used in function in which analyses are
#' done for all matches between these teams. This is done by loading the saved .RData for
#' each match and performing an rbind of the data frames for each match
#'
#' @usage
#' getAllMatchesBetweenTeams(team1,team2,dir=".",save=FALSE)
#'
#' @param team1
#' One of the team in consideration e.g (India, Australia, England)
#'
#' @param team2
#' The other team for which matches are needed e.g( India, Sri Lanka, Pakistan)
#'
#' @param dir
#' The directory which has the RData files of matches between teams
#'
#' @param save
#' Default=FALSE. This parameter indicates whether the combined data frame needs to be saved or not. It is recommended
#' to save this large dataframe as the creation of this data frame takes a several seconds depending
#' on the number of matches
#'
#' @return matches
#' The combined data frame
#'
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' # Get all matches for team India
#' getAllMatchesAllOpposition("India",dir="../data/",save=TRUE)
#' getAllMatchesAllOpposition("Australia",dir="./mysavedata/",save=TRUE)
#' }
#'
#' @seealso
#' \code{\link{bowlerMovingAverage}}\cr
#' \code{\link{bowlerWicketPlot}}\cr
#' \code{\link{bowlerWicketsVenue}}\cr
#' \code{\link{getAllMatchesAllOpposition}}\cr
#'
#' @export
#'

getAllMatchesBetweenTeams <- function(team1,team2,dir=".",save=FALSE){
     overs=NULL
    # Create 2 filenames with both combinations of team1 and team2
    d1 <- paste(team1,"-",team2,"*",sep="")
    d2 <- paste(team2,"-",team1,"*",sep="")
    path1=paste(dir,"/",d1,sep="")
    path2=paste(dir,"/",d2,sep="")
    # Capture both combinations
    fl1 <- Sys.glob(path1)
    fl2 <- Sys.glob(path2)
    fl3 <-c(fl1,fl2)
    if(length(fl3) != 0){
    
    # Create a data frame with all matches
    matches <- NULL
    for(i in 1:length(fl3)){
        load(fl3[i])
        matches <- rbind(matches,overs)
    }
    b <- paste(team1,"-",team2,"-allMatches.RData",sep="")
    if(save){
        save(matches,file=b)
    }

    matches
    }
}
