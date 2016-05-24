##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 20 Mar 2016
# Function: convertYaml2RDataframe
# This function converts a given yaml file to a data frame and stores this as n .RData
# The yaml file is read from a given source directory coverted and saved to a target directory.
# The target file is the created by using the name of the opposing teams and the date of the match
#
###########################################################################################
#' @title
#' Converts and save yaml files to dataframes
#'
#' @description
#' This function coverts all Yaml files from source directory to data frames. The data frames
#' are then stored as .RData. The saved file is of the format team1-team2-date.RData
#' For e.g. England-India-2008-04-06.RData etc
#' @usage
#' convertYaml2RDataframe(yamlFile,sourceDir=".",targetDir=".")
#'
#' @param yamlFile
#' The  yaml file to be converted to dataframe and saved
#'
#' @param sourceDir
#' The source directory of the yaml file
#'
#' @param targetDir
#' The target directory in which the data frame is stored as RData file
#'
#' @return None
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
#' # In the example below ../yamldir c
#' convertYaml2RDataframe("225171.yaml",".","../data")
#' }
#'
#' @seealso
#' \code{\link{bowlerMovingAverage}}\cr
#' \code{\link{bowlerWicketPlot}}\cr
#' \code{\link{bowlerWicketsVenue}}\cr
#' \code{\link{convertAllYaml2RDataframes}}\cr
#'
#' @export
#'
convertYaml2RDataframe <- function(yamlFile,sourceDir=".",targetDir="."){

    yaml.load_file=info.dates=info.match_type=info.overs=info.venue=NULL
    info.teams=matchType=winner=result=venue=NULL
    pth = paste(sourceDir,"/",yamlFile,sep="")
    print(pth)

    # Load yaml file
    a <- yaml.load_file(pth)

    # Cast as data frame for easy processing
    b <- as.data.frame(a)
    sz <- dim(b)

    # Gather the meta information
    meta <- select(b,info.dates,info.match_type,info.overs, info.venue,
                   info.teams)
    names(meta) <- c("date","matchType","overs","venue","team1")

    # Check if there was a winner or if there was no result (tie,draw)
    if(!is.null(b$info.outcome.winner)){
        meta$winner <- b$info.outcome.winner
        meta$result <- "NA"
    } else if(!is.null(b$info.result)){
        meta$winner <- "NA"
        meta$result <- b$info.result
    } else if(!is.null(b$info.outcome.result)){
        meta$winner <- "NA"
        meta$result <- b$info.outcome.result
    }

    meta$team2 = meta[2,5]
    meta <- meta[1,]


    #Reorder columns
    meta <- select(meta,date,matchType,overs,venue,team1,team2,winner,result)

    # Remove the innings and deliveries from the column names
    names(b) <-gsub("innings.","",names(b))
    names(b) <- gsub("deliveries.","",names(b))

    # Create an empty data frame
    overs <- data.frame(ball=character(),team=factor(),batsman=factor(),
                        bowler=factor(),nonStriker=factor(),byes=numeric(),
                        legbyes=numeric(), noballs=numeric(), wides=numeric(),
                        nonBoundary=factor(), penalty=factor(),
                        runs=factor(),extras=factor(),totalRuns=factor(),
                        wicketFielder=character(), wicketKind=character(),
                        wicketPlayerOut=character(),date=factor(),
                        matchType=factor(),
                        overs=integer(),venue=factor(),team1=factor(),team2=factor(),
                        winner=character(),result=character())


    #Choose the columns which have the ball by ball detail
    idx = which(names(b) == "1st.0.1.batsman")
    m <- b[1,idx:sz[2]]
    #Transpose to the details of the match
    match <- t(m)
    rnames <- rownames(match)

    match <- as.data.frame(cbind(rnames,match))

    # Gather details for first team

    # Set the number of overs to 50 for ODI matches
    numOver <- seq(from=0,to=50,by=1)

    # Create string of delivery in each over upto delivery 16 in case of no balls,
    # wides etc.
    # Note: The over can be more than .6 when you have no balls, wides etc
    d <- c(".1",".2",".3",".4",".5",".6",".7",".8",".9",".11",".12",
           ".13",".14",".15",".16")

    m <- 1
    # Create a vector of deliveries from 0 to 50 by concatenating string
    delivery <- NULL
    for(k in 1:length(numOver)){
        for(l in 1:length(d)){
            delivery[m] <- paste(numOver[k],d[l],sep="")
            m=m+1
        }
    }


    #Create string for 1st team
    print("first loop")
    s <- paste("1st.",delivery,"\\.",sep="")
    team1 <- b$`1st.team`[1]
    # Parse the yaml file over by over and store as a row of data
    overs1 <- parseYamlOver(match,s,team1,overs,delivery,meta)


    # Create string for 2nd team
    print("second loop")
    s1 <- paste("2nd.",delivery,"\\.",sep="")
    team2 <- b$`2nd.team`[1]
    overs2 <- parseYamlOver(match,s1,team2,overs,delivery,meta)

    # Row bind the 1dst
    overs <- rbind(overs1,overs2)

    # Change factors to appropiate type
    overs$byes <- as.numeric(as.character(overs$byes))
    overs$legbyes <- as.numeric(as.character(overs$legbyes))
    overs$wides <- as.numeric(as.character(overs$wides))
    overs$noballs <- as.numeric(as.character(overs$noballs))
    overs$nonBoundary <- as.numeric(as.character(overs$nonBoundary))
    overs$penalty <- as.numeric(as.character(overs$penalty))
    overs$runs <- as.numeric(as.character(overs$runs))
    overs$extras <- as.numeric(as.character(overs$extras))
    overs$totalRuns <- as.numeric(as.character(overs$totalRuns))
    overs$date = as.Date(overs$date)
    overs$overs <- as.numeric(as.character(overs$overs))
    sapply(overs,class)


    teams <- as.character(unique(overs$team))
    #Create a unique file which is based on the opposing teams and the date of the match
    filename <- paste(meta$team1,"-",meta$team2,"-",meta$date,".",
                      "RData",sep="")

    to <- paste(targetDir,"/",filename,sep="")
    # Save as .RData
    save(overs,file=to)
    # Write the name of the file that was converted and the converted file for reference
    convertedFile <- paste(yamlFile,filename,sep=":")
    outputFile <- paste(targetDir,"/","convertedFiles.txt",sep="")
    write(convertedFile,outputFile,append=TRUE)


}




