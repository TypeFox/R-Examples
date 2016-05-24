###################################################
### chunk number 1: 
###################################################



# 'readData' to reads the data of a specified disease of several years
#       and generates a state chain using the bulletin knowledge
#
# Parameter:
#       abb : abbreviation of the disease
#       week53to52: Boolean indicating whether to convert RKI 53 Weeks System to 52 weeks a year
readData <- function(abb,week53to52=TRUE,sysPath=TRUE){
  #Read depending on which path is requested
  if (sysPath) {
    #Prepend the systempath/data to the filename
    #hoehle 2012-07-24 - this does not work when package is not
    #installed. Use extdata as recommended in the file package structure.
    file <- file.path(path.package('surveillance'),'extdata',paste(abb,".txt",sep=""))
  } else {
    file <- file.path(paste(abb,".txt",sep=""))
  }

  # read the data from four years and write it to a table
  #file <- paste( dataPath, abb , ".txt" , sep="" )
  fileTable <- read.table( file=file, header=TRUE )
  observed <- fileTable$observed
  state <- fileTable$state
  
  result = list(observed=observed, state=state)
  
  class(result) = "disProg" # for disease progress
  
  #Convert to 52 week system...
  if (week53to52) {
    result <- correct53to52(result)
  }
  
  result$freq <- 52
  result$start <- c(2001,1)
  
  return(result)
}




###################################################
### chunk number 2: 
###################################################

toFileDisProg <- function(disProgObj, toFile){

        length <- length(disProgObj$observed)

        writeMatrix <- matrix(0, length, 3)
        dimnames(writeMatrix) <- list(c(), c("week", "observed", "state"))

        writeMatrix[,"week"] <- 1:length
        writeMatrix[,"observed"] <- disProgObj$observed
        writeMatrix[,"state"] <- disProgObj$state

        write.table(writeMatrix, toFile, row.names = FALSE, sep = "\t")
}




###################################################
### chunk number 3: 
###################################################

# 'correct53to52' sums up and cuts a value from a splited last and first week of a year
#
# Parameter:
#       disProgObj - object of class disProgObj (including the observed and the state chain)
#       firstweek: the number in observed of the first week in a year, default = 1
# ouput:
#       disProgObj: the new disProgObj


correct53to52 <- function(disProgObj, firstweek = 1){

        if(firstweek > length(disProgObj$observed)){
                stop("firstweek doesn't exist")
        }

        observed <- disProgObj$observed
        state <- disProgObj$state

        if(length(state) != length(observed)){
                stop("state and observed don't have the same length")
        }

        # do not cut, if observed is too short
        length = length(observed[firstweek:length(observed)])

        if(length > 53){

                lastyear <- floor((length-1)/53)
                # sum case numbers of double weeks up
                for(i in 1:lastyear){
                        # last week of year i (-i+1 because the array now is shorter)
                        last <- firstweek + i * 52
                        # first week in year i+1
                        firstnew <- last + 1
                        observed[firstnew]  <- observed[last]  + observed[firstnew]
                        # delete double weeks
                        observed <- observed[-c(last)]

                        # with state
                        state[firstnew]  <- state[last]  + state[firstnew]
                        # delete double weeks
                        state <- state[-c(last)]
                }
        }

        # correct also the first week, if it doesn't is the beginning
        if(firstweek > 1){
                observed[firstweek] <- observed[firstweek] + observed[firstweek-1]
                observed <- observed[-c(firstweek-1)]
                state[firstweek] <- state[firstweek] + state[firstweek-1]
                state <- state[-c(firstweek-1)]
        }

        # correct all 2 to 1
        state[state==2] <- 1

        disProgObj$observed <- observed
        disProgObj$state <- state

        return(disProgObj)
}


###################################################
### chunk number 4: 
###################################################

enlargeData <- function(disProgObj, range = 1:156, times = 1){

        # enlarge observed
        disProgObj$observed <- c(rep(disProgObj$observed[range], times), disProgObj$observed)
        # enlarge state
        disProgObj$state <- c(rep(disProgObj$state[range], times), disProgObj$state)

        return(disProgObj)
}



