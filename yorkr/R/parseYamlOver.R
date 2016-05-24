##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 20 Mar 2016
# Function: parseYamlOver
# This function converts a given yaml file to a dataframe delivery by delivery.
#
###########################################################################################
#' @title
#' Parse yaml file and convert to dataframe
#'
#' @description
#' This function parses the yaml file and converts it into a data frame. This is an internal function and
#' is used by convertAllYaml2RDataframes() & convertYaml2RDataframe()
#'
#' @usage
#' parseYamlOver(match,s,ateam,overs,delivery,meta)
#'
#' @param match
#' The dataframe of the match
#'
#' @param s
#' The string with the delivery
#'
#' @param ateam
#' The team
#'
#' @param overs
#' overs
#'
#' @param delivery
#' The delivery of the over
#'
#' @param meta
#' The meta information of the match
#'
#' @return overs
#' The dataframe of overs
#'
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#'
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
#'
parseYamlOver <- function(match,s,ateam,overs,delivery,meta) {
    team=ball=totalRuns=rnames=batsman=bowler=nonStriker=NULL
    byes=legbyes=noballs=wides=nonBoundary=penalty=runs=NULL
    extras=wicketFielder=wicketKind=wicketPlayerOut=NULL
    #The xsflag will be used when there are more than 10 deliveries in the over
    xsflag=FALSE

    # Loop through all deliveries one by one. Check whether the over had 6,7,8,9,10 or more deliveries
    for(i in 1:length(delivery)){

        # Filter rows based on the delivery(ball) as overset
        overset <- filter(match,grepl(s[i],rnames))
        #Transpose
        over <- as.data.frame(t(overset))
        # Generate a row vector
        over <- over[2,]

        #Check the number of deliveries in the over
        d <- dim(over)


        # The over had 6 deliveries
        # Check the number of deliveries in the over
        # If the number of deliveries in the over was
        # 6 : This is a notmal over. No extras, no wickets
        # 7 : There was 1 wide/no ball/legbye/bye/non-boundary/penalty
        # 8 : There was a wickt in the over (possibly bowled)
        # 9 : There was a wicket in this over. wicket, wicketKind and wicketPlayerOut
        # 10 : There was a bye/legbyes/wide/noball etc and a player got out
        # 10+ : There was more than 10 balls. Split into 10 + x

        if(d[2] == 6){
            # Set the names of the columns
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
            over$ball=gsub("\\\\.","",s[i])
            over$team = ateam

            # Reorder the  columns
            over <- select(over, ball,team,batsman,bowler,nonStriker,
                           byes,legbyes,noballs,
                           wides,nonBoundary,penalty, runs,
                           extras,totalRuns,wicketFielder,
                           wicketKind,wicketPlayerOut)

            over <- cbind(over,meta)

        } else if(d[2]==7){
            # The over had 7 deliveries. The extra delivery is because of a
            # legbye, bye, wide, no ball,non-boundary or penalty.
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
                cat("parse1=",i,"\n")
                names(over) <-c("batsman","bowler","nonStriker","runs","extras","nonBoundary","totalRuns")
                over$byes=as.factor(0)
                over$legbyes=as.factor(0)
                over$wides=as.factor(0)
                over$noballs=as.factor(0)
                over$penalty=as.factor(0)
            } else if(sum(grepl("penalty",overset$rnames))){
                names(over) <-c("batsman","bowler","penalty","nonStriker","runs","extras","totalRuns")
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
            over$ball=gsub("\\\\.","",s[i])
            over$team = ateam
            # Reorder
            over <- select(over, ball,team,batsman,bowler, nonStriker,
                           byes,legbyes,noballs,
                           wides,nonBoundary,penalty, runs,
                           extras,totalRuns,wicketFielder,
                           wicketKind,wicketPlayerOut)
            #over <- over[,c(14,15,1,2,3,4,5,6,7,8,9,10,11,12,13)]
            over <- cbind(over,meta)
        } else if(d[2] ==8){ # A player got out in this over
            names(over) <-c("batsman","bowler","nonStriker","runs","extras","totalRuns",
                            "wicketKind","wicketPlayerOut")

            # Add the missing elements for extras
            over$byes<-as.factor(0)
            over$legbyes<-as.factor(0)
            over$noballs<-as.factor(0)
            over$wides<-as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty<-as.factor(0)

            over$wicketFielder="nobody"
            over$ball=gsub("\\\\.","",s[i])
            over$team = ateam

            # Reorder the columns
            over <- select(over, ball,team,batsman,bowler,nonStriker,
                           byes,legbyes,noballs,
                           wides,nonBoundary,penalty,runs,
                           extras,totalRuns,wicketFielder,
                           wicketKind,wicketPlayerOut)
            #over <- over[,c(14,15,1,2,9,10,11,12,3,4,5,6,10,7,8)]
            over <- cbind(over,meta)


        } else if(d[2] ==9){ # A player got out in this over
            names(over) <-c("batsman","bowler","nonStriker","runs","extras","totalRuns",
                            "wicketFielder","wicketKind","wicketPlayerOut")

            # Add the missing elements for extras
            over$byes<-as.factor(0)
            over$legbyes<-as.factor(0)
            over$noballs<-as.factor(0)
            over$wides<-as.factor(0)
            over$nonBoundary <- as.factor(0)
            over$penalty<-as.factor(0)

            over$ball=gsub("\\\\.","",s[i])
            over$team = ateam
            over <- select(over, ball,team,batsman,bowler,nonStriker,
                           byes,legbyes,noballs,
                           wides,nonBoundary,penalty, runs,
                           extras,totalRuns,wicketFielder,
                           wicketKind,wicketPlayerOut)
            over <- cbind(over,meta)

        } else if(d[2] == 10) { # The a player got out in this over
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

                print("OOOOOO*************************************************")
                break
                over$byes=as.factor(0)
                over$legbyes=as.factor(0)
                over$wides=as.factor(0)
                over$noballs=as.factor(0)
                over$penalty=as.factor(0)
            } else if(sum(grepl("penalty",overset$rnames))){
                names(over) <-c("batsman","bowler","penalty","nonStriker","runs","extras","totalRuns")
                over$byes=as.factor(0)
                over$legbyes=as.factor(0)
                over$noballs=as.factor(0)
                over$wides=as.factor(0)
                over$nonBoundary <- as.factor(0)
            }

            over$ball=gsub("\\\\.","",s[i])
            over$team = ateam

            # Reorder the columns
            over <- select(over, ball,team,batsman,bowler,nonStriker,
                           byes,legbyes,noballs,
                           wides,nonBoundary,penalty,runs,
                           extras,totalRuns,wicketFielder,
                           wicketKind,wicketPlayerOut)
            over <- cbind(over,meta)

        }else if(d[2] == 0){
            next
        } else if(d[2] >10){
            # This situation can arise when the the number of no balls and
            # wides take it over the 10 deliveries.
            # We need to split the deliveries intp 2 groups - group1 & group2
            # The two group are required because the deliveries
            # start from 0.1,0.2,0.3 ...0.9. There is no 0.10, instead this is
            # written as over.delivery.*.1. The usual nomeclature is
            # over.delivery.batsman(bowler,etc) e.g 2nd.2.1.batsman,
            # 2nd.2.1.bowler and so on. The 10th delivery is written as
            # 2nd.2.1.batsman.1,2nd.2nd.1.bowler.1
            # Group1
            overset <- filter(match,grepl(s[i],rnames))
            jj<-paste(s[i],"batsman",sep="")
            ll <- grepl(jj,overset$rnames)
            mm <- which(ll)
            dist1 = mm[2] - mm[1]
            o1 <- overset[1:dist1,]
            o1 <- as.data.frame(t(o1))
            o1 <- o1[2,]
            d <- dim(o1)
            overinit <- specialProc(d[2],overset,ateam,o1,s[i],meta)

            # Group 2
            o2 <- overset[mm[2]:length(ll),]
            o2 <- as.data.frame(t(o2))
            o2 <- o2[2,]
            d <- dim(o2)
            overfinal <- specialProc(d[2],overset,ateam,o2,s[i],meta)

            xsflag=TRUE

        }

        # Row bind the data
        #print(dim(overs))
        #print(dim(over))

        if(!xsflag){
            overs <- rbind(overs,over)
        }

        if(xsflag){
            overs <- rbind(overs,overinit)
            overs <- rbind(overs,overfinal)
            xsflag=FALSE
        }
    }
    overs
}
