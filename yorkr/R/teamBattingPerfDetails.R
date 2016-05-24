##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 25 Mar 2016
# Function: teamBattingPerfDetails
# This function get the overall team batting details of the matcjh
#
###########################################################################################
#' @title
#' Gets the team batting details.
#'
#' @description
#' This function gets the team batting detals
#'
#' @usage
#' teamBattingPerfDetails(match,theTeam,includeInfo=FALSE)
#'
#' @param match
#' The match between the teams
#'
#' @param theTeam
#' The team for which the the batting partnerships are sought
#'
#' @param includeInfo
#' Whether to include venue,date, winner and result
#'
#' @return df
#' dataframe
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
#' #teamBattingPerfDetails()
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBattingScorecardMatch}}\cr
#'
#'
#'
teamBattingPerfDetails <- function(match,theTeam,includeInfo=FALSE){
    team=batsman=runs=fours=sixes=NULL
    byes=legbyes=noballs=wides=bowler=wicketFielder=NULL
    wicketKind=wicketPlayerOut=NULL
    # Initialise to NULL
    details <- NULL
    a <-filter(match,team==theTeam)
    sz <- dim(a)
    if(sz[1] == 0){
        #cat("No batting records.\n")
        return(NULL)
    }
    b <- select(a,batsman,runs)
    names(b) <-c("batsman","runs")

    #Compute the number of 4s
    c <-
        b %>%
        mutate(fours=(runs>=4 & runs <6)) %>%
        filter(fours==TRUE)

    # Group by batsman. Count 4s
    d <-    summarise(group_by(c, batsman),fours=n())

    # Get the total runs for each batsman
    e <-summarise(group_by(a,batsman),sum(runs))
    names(b) <-c("batsman","runs")
    details <- full_join(e,d,by="batsman")
    names(details) <-c("batsman","runs","fours")

    # Compute number of 6's
    f <-
        b %>%
        mutate(sixes=(runs ==6)) %>%
        filter(sixes == TRUE)
    # Group by batsman. COunt 6s
    g <- summarise(group_by(f, batsman),sixes=n())
    names(g) <-c("batsman","sixes")

    # Full join with 4s and 6s
    details <- full_join(details,g,by="batsman")

    # Count the balls played by the batsman
    ballsPlayed <-
        a  %>%
        select(batsman,byes,legbyes,wides,noballs,runs) %>%
        filter(wides ==0,noballs ==0,byes ==0,legbyes == 0) %>%
        select(batsman,runs)

    ballsPlayed<- summarise(group_by(ballsPlayed,batsman),count=n())
    names(ballsPlayed) <- c("batsman","ballsPlayed")

    # Create a data frame
    details <- full_join(details,ballsPlayed,by="batsman")

    # If there are NAs then replace with 0's
    if(sum(is.na(details$fours)) != 0){
        details[is.na(details$fours),]$fours <- 0
    }
    if(sum(is.na(details$sixes)) != 0){
        details[is.na(details$sixes),]$sixes <- 0
    }

    details <- select(details,batsman,ballsPlayed,fours,sixes,runs)



    #Calculate strike rate
    details <- mutate(details,strikeRate=round(((runs/ballsPlayed)*100),2))

    w <- filter(a,wicketKind !="not-out" | wicketPlayerOut != "nobody" )

    # Remove unnecessary factors
    w$wicketPlayerOut <-factor(w$wicketPlayerOut)

    wkts <- select(w,batsman,bowler,wicketFielder,wicketKind,wicketPlayerOut)

    details <- full_join(details,wkts,by="batsman")

    # Set as character to be able to assign value
    details$wicketPlayerOut <- as.character(details$wicketPlayerOut)
    details$wicketKind <- as.character(details$wicketKind)
    details$wicketFielder <- as.character(details$wicketFielder)
    details$bowler <- as.character(details$wicketFielder)

    # Set the NA columns in wicketPlayerOut with notOut
    # Also set the other columns for this row
    if(sum(is.na(details$wicketPlayerOut))!= 0){
        details[is.na(details$wicketPlayerOut),]$wicketPlayerOut="notOut"
        details[is.na(details$wicketKind),]$wicketKind="notOut"
        details[is.na(details$wicketFielder),]$wicketFielder="nobody"
        details[is.na(details$bowler),]$bowler="nobody"
    }

    # Determine the opposition
    t <- match$team != theTeam
    # Pick the 1st element

    t1 <- match$team[t]
    opposition <- as.character(t1[1])

    if(includeInfo == TRUE) {
        details$date <- a$date[1]
        details$venue <- a$venue[1]
        details$opposition <- opposition
        details$winner <- a$winner[1]
        details$result <- a$result[1]
    }


    details

}
