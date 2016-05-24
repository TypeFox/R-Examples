#------------------------------------------------------------------------------
#' modifiedRecursive trimming procedure.
#'
#' \code{modifiedRecursive} takes a data frame of RT data and returns trimmed rt
#' data that fall below a set standard deviation above the each participant's
#' mean for each condition, with the criterion changing as more trials are
#' removed,  as described in van Selst & Jolicoeur (1994).
#' @param data A data frame. It must contain columns named "participant",
#' "condition", "rt", and "accuracy". The RT can be in seconds
#' (e.g., 0.654) or milliseconds (e.g., 654). Typically, "condition" will
#' consist of strings. "accuracy" must be 1 for correct and 0 for error
#' responses.
#' @param minRT The lower criteria for acceptable response time. Must be in
#' the same form as rt column in data frame (e.g., in seconds OR milliseconds).
#' All RTs below this value are removed before proceeding with SD trimming.
#' @param omitErrors If set to TRUE, error trials will be removed before
#' conducting trimming procedure. Final data returned will not be influenced
#' by errors in this case.
#' @param digits How many decimal places to round to after trimming?
#'
#' @references Van Selst, M. & Jolicoeur, P. (1994). A solution to the effect
#' of sample size on outlier elimination. \emph{Quarterly Journal of Experimental
#' Psychology, 47} (A), 631-650.
#'
#'
#' @examples
#' # load the example data that ships with trimr
#' data(exampleData)
#'
#' # perform the trimming, returning mean RT
#' trimmedData <- modifiedRecursive(data = exampleData, minRT = 150)
#'
#' @importFrom stats sd
#'
#' @export

modifiedRecursive <- function(data, minRT, omitErrors = TRUE, digits = 3){


  # remove errors if the user has asked for it
  if(omitErrors == TRUE){
    trimmedData <- subset(data, data$accuracy == 1)
  } else {
    trimmedData <- data
  }

  # get the list of participant numbers
  participant <- sort(unique(trimmedData$participant))

  # get the list of experimental conditions
  conditionList <- unique(trimmedData$condition)

  # trim the data to remove trials below minRT
  trimmedData <- subset(trimmedData, trimmedData$rt > minRT)

  # ready the final data set
  finalData <- matrix(0, nrow = length(participant),
                      ncol = length(conditionList))

  # give the columns the condition names
  colnames(finalData) <- conditionList

  # add the participant column
  finalData <- cbind(participant, finalData)

  # convert to data frame
  finalData <- data.frame(finalData)

  # intialise looping variable for subjects
  i <- 1

  # loop over all subjects
  for(currSub in participant){

    # intialise looping variable for conditions. It starts at 2 because the
    # first column in the data file containing condition information is the
    # second one.
    j <- 2

    # loop over all conditions
    for(currCond in conditionList){

      # get the relevant data
      tempData <- subset(trimmedData, trimmedData$participant == currSub &
                           trimmedData$condition == currCond)


      # find the average, and add to the data frame
      finalData[i, j] <- round(modifiedRecursiveTrim(tempData$rt),
                               digits = digits)

      # update condition loop counter
      j <- j + 1
    }

    # update participant loop counter
    i <- i + 1
  }
  return(finalData)
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### The function that acually does the trimming
modifiedRecursiveTrim <- function(data){

  # load the linear interpolation data file (hidden from user)
  criterion <- linearInterpolation

  # repeat the following
  repeat{


    # the data needs to have more than 2 RTs in order to function.
    # (thanks to Ayala Allon for noticing this.)
    if(length(data) <= 2){
      break
    }

    ### get the parameters for the moving criterion

    # get the sample size of the data
    sampleSize <- length(data)

    # if the sample size is greater than 100, use SDs for N = 100
    if(sampleSize > 100){
      sampleSize <- 100
    }

    # look up the SD to use for the current sampleSize
    stDev <- criterion$nonRecursive[sampleSize]

    ### now do the removal of trials

    # what is the maximum RT value in the current data?
    x <- max(data)

    # temporarily remove it
    tempData <- data[data != x]

    # find the mean, SD, & cutoffs for temporary data
    curMean <- mean(tempData)
    curSD <- sd(tempData)
    maxCutoff <- curMean + (stDev * curSD)
    minCutoff <- curMean - (stDev * curSD)

    # now, go back to the main data (with temporary RTs replaced)

    # find the largest & smallest RT values
    x <- max(data)
    y <- min(data)

    # initalise this variable before looping. This will store how many trials
    # have been removed on each recursion. The loop stops when no trials
    # have been removed on a particular iteration
    removedTrials <- 0

    # if there is a data point above the cutoff, then remove it
    if(x > maxCutoff){
      data <- data[data != x]
      removedTrials <- 1
    }

    # if there is a data point below the cutoff, then remove it
    if(y < minCutoff){
      data <- data[data != y]
      removedTrials <- 1
    }

    # when there are no trials removed on the current iteration, break the loop
    if(removedTrials == 0){
      break
    }

    # alternatively, stop the loop when the data set hits 4 RTs left (as in
    # vanSelst & Jolicoeur, 1994)
    if(length(data) < 5){
      break
    }

  } # end of repeat loop

  # now compute the mean of the final data set. If there are only two RTs,
  # then return NA
  if(length(data) == 2){
    finalData <- NA
  } else {
    finalData <- mean(data)
  }

  # return the data
  return(finalData)

} # end of function

#------------------------------------------------------------------------------
