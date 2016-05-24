##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 20 Mar 2016
# Function: specialProc
# This is a helper function used by parseYamlOver when the over has more tan 10
# deliveries.
#
###########################################################################################
#' @title
#' Used to parse yaml file
#'
#' @description
#' This is special processing function. This is an internal function and
#' is used by convertAllYaml2RDataframes() & convertYaml2RDataframe()
#'
#' @usage
#' specialProc(dist, overset, ateam,over,str1,meta)
#'
#' @param dist
#' dist
#'
#' @param overset
#' overset
#'
#' @param ateam
#' The team
#'
#' @param over
#' over
#'
#' @param str1
#' str1
#'
#'
#' @param meta
#' The meta information of the match
#'
#' @return over
#' The dataframe of over
#'
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#'
#' @author
#' Tinniam V Ganesh
#'
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' # Parse the yaml over
#' }
#'
#' @seealso
#' \code{\link{getBatsmanDetails}}\cr
#' \code{\link{getBowlerWicketDetails}}\cr
#' \code{\link{batsmanDismissals}}\cr
#' \code{\link{getTeamBattingDetails}}\cr
#'
#'
# This functio is used when there are more than 10 deliveries in the over
specialProc <- function(dist, overset, ateam,over,str1,meta){
    team=ball=totalRuns=rnames=batsman=bowler=nonStriker=i=NULL
    byes=legbyes=noballs=wides=nonBoundary=penalty=runs=NULL
    extras=wicketFielder=wicketKind=wicketPlayerOut=NULL
    if(dist == 6){
        names(over) <-c("batsman","bowler","nonStriker","runs","extras","totalRuns")

        # Add the missing elements for extras
        over$byes<-as.factor(0)
        over$legbyes<-as.factor(0)
        over$noballs<-as.factor(0)
        over$wides<-as.factor(0)
        over$nonBoundary <- as.factor(0)
        over$penalty<-as.factor(0)

        over$wicketFielder="nobody"
        over$wicketKind="not-out"
        over$wicketPlayerOut="nobody"
        over$ball=gsub("\\\\.","",str1)
        over$team = ateam
        # Reorder the rows
        over <- select(over, ball,team,batsman,bowler,nonStriker,
                       byes,legbyes,noballs,
                       wides,nonBoundary,penalty,runs,
                       extras,totalRuns,wicketFielder,
                       wicketKind,wicketPlayerOut)

        over <- cbind(over,meta)

    } else if(dist==7){
        # The over had 7 deliveries
        if(sum(grepl("\\.byes",overset$rnames))){
            names(over) <-c("batsman","bowler","byes","nonStriker","runs","extras","totalRuns")
            over$legbyes=as.factor(0)
            over$noballs=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty=as.factor(0)
        } else if(sum(grepl("legbyes",overset$rnames))){
            names(over) <-c("batsman","bowler","legbyes","nonStriker","runs","extras","totalRuns")
            over$byes=as.factor(0)
            over$noballs=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty=as.factor(0)
        } else if(sum(grepl("noballs",overset$rnames))){
            names(over) <-c("batsman","bowler","noballs","nonStriker","runs","extras","totalRuns")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty=as.factor(0)
        } else if(sum(grepl("wides",overset$rnames))){
            names(over) <-c("batsman","bowler","wides","nonStriker","runs","extras","totalRuns")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$noballs=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty=as.factor(0)
        } else if(sum(grepl("non_boundary",overset$rnames))){
            cat("sp=",i,"\n")
            names(over) <-c("batsman","bowler","nonStriker","runs","extras","nonBoundary","totalRuns")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$wides=as.factor(0)
            over$noballs=as.factor(0)
            over$penalty=as.factor(0)
        } else if(sum(grepl("penalty",overset$rnames))){
            names(over1) <-c("batsman","bowler","penalty","nonStriker","runs","extras","totalRuns")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$noballs=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
        }

        # Add missing elements
        over$wicketFielder="nobody"
        over$wicketKind="not-out"
        over$wicketPlayerOut="nobody"
        over$ball=gsub("\\\\.","",str1)
        over$team = ateam
        # Reorder
        over <- select(over, ball,team,batsman,bowler,nonStriker,
                       byes,legbyes,noballs,
                       wides,nonBoundary,penalty,runs,
                       extras,totalRuns,wicketFielder,
                       wicketKind,wicketPlayerOut)
        #over <- over[,c(14,15,1,2,3,4,5,6,7,8,9,10,11,12,13)]
        over <- cbind(over,meta)
        #cat("Hhhh",dim(over),"\n")
    } else if(dist ==8){
        names(over) <-c("batsman","bowler","nonStriker","runs","extras","totalRuns","wicketKind","wicketPlayerOut")

        # Add the missing elements for extras
        over$byes<-as.factor(0)
        over$legbyes<-as.factor(0)
        over$noballs<-as.factor(0)
        over$wides<-as.factor(0)
        over$nonBoundary <- as.factor(0)
        over$penalty<-as.factor(0)

        over$wicketFielder="nobody"
        over$ball=gsub("\\\\.","",str1)
        over$team = ateam
        # Reorder
        over <- select(over, ball,team,batsman,bowler,nonStriker,
                       byes,legbyes,noballs,
                       wides,nonBoundary,penalty,runs,
                       extras,totalRuns,wicketFielder,
                       wicketKind,wicketPlayerOut)
        over <- cbind(over,meta)


    } else if(dist ==9){
        names(over) <-c("batsman","bowler","nonStriker","runs","extras","totalRuns",
                        "wicketFielder","wicketKind","wicketPlayerOut")

        # Add the missing elements for extras
        over$byes<-as.factor(0)
        over$legbyes<-as.factor(0)
        over$noballs<-as.factor(0)
        over$wides<-as.factor(0)
        over$nonBoundary <- as.factor(0)
        over$penalty<-as.factor(0)

        over$ball=gsub("\\\\.","",str1)
        over$team = ateam
        over <- select(over, ball,team,batsman,bowler,nonStriker,
                       byes,legbyes,noballs,
                       wides,nonBoundary,penalty,runs,
                       extras,totalRuns,wicketFielder,
                       wicketKind,wicketPlayerOut)
        over <- cbind(over,meta)


    } else if(dist == 10){
        if(sum(grepl("\\.byes",overset$rnames))){
            names(over) <-c("batsman","bowler","byes","nonStriker","runs","extras","totalRuns",
                            "wicketFielder","wicketKind","wicketPlayerOut")
            over$legbyes=as.factor(0)
            over$noballs=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty<-as.factor(0)
        } else if(sum(grepl("legbyes",overset$rnames))){
            names(over) <-c("batsman","bowler","legbyes","nonStriker","runs","extras","totalRuns",
                            "wicketFielder","wicketKind","wicketPlayerOut")
            over$byes=as.factor(0)
            over$noballs=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty<-as.factor(0)
        } else if(sum(grepl("noballs",overset$rnames))){
            names(over) <-c("batsman","bowler","noballs","nonStriker","runs","extras","totalRuns",
                            "wicketFielder","wicketKind","wicketPlayerOut")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty<-as.factor(0)
        } else if(sum(grepl("wides",overset$rnames))){
            names(over) <-c("batsman","bowler","wides","nonStriker","runs","extras","totalRuns",
                            "wicketFielder","wicketKind","wicketPlayerOut")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$noballs=as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty<-as.factor(0)
        } else if(sum(grepl("non_boundary",overset$rnames))){
            cat("sp=",i,"\n")
            names(over) <-c("batsman","bowler","nonStriker","runs","extras","nonBoundary","totalRuns")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$wides=as.factor(0)
            over$noballs=as.factor(0)

            over$penalty=as.factor(0)
            over$penalty<-as.factor(0)
        } else if(sum(grepl("penalty",overset$rnames))){
            names(over) <-c("batsman","bowler","penalty","nonStriker","runs","extras","totalRuns")
            over$byes=as.factor(0)
            over$legbyes=as.factor(0)
            over$noballs=as.factor(0)
            over$wides=as.factor(0)
            over$nonBoundary=as.factor(0)
            over$penalty<-as.factor(0)
        }

        over$ball=gsub("\\\\.","",str1)
        over$team = ateam
        over <- select(over, ball,team,batsman,bowler,nonStriker,
                       byes,legbyes,nonBoundary,penalty,noballs,
                       wides,runs,
                       extras,totalRuns,wicketFielder,
                       wicketKind,wicketPlayerOut)

        over <- cbind(over,meta)
        print("Ho!")
    }
    #cat("returning",dim(over),"\n")
    over

}
