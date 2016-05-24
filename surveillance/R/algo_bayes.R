###################################################
### chunk number 1:
###################################################


# Implementation of the Bayes system.
# The system evaluates specified timepoints and gives alarm if it recognizes
# an outbreak for this timepoint.
#
# Features:
# Choice between different Bayes sub-systems (difference in reference values).

algo.bayesLatestTimepoint <- function(disProgObj, timePoint = NULL, control = list(b = 0, w = 6, actY = TRUE, alpha=0.05)){

  observed <- disProgObj$observed
  freq <- disProgObj$freq

  # If there is no value in timePoint, then take the last value in observed
  if(is.null(timePoint)){
        timePoint = length(observed)
  }

  #If no level specified.

  # check if the vector observed includes all necessary data.
  if((timePoint-(control$b*freq)-control$w) < 1){
        stop("The vector of observed is too short!")
  }

  # construct the reference values
  basevec <- c()
  # if actY == TRUE use also the values of the year of timepoint
  if(control$actY){
        basevec <- observed[(timePoint - control$w):(timePoint - 1)]
  }
  # check if you need more referencevalues of the past
  if(control$b >= 1){
    for(i in 1:control$b){
        basevec <- c(basevec, observed[(timePoint-(i*freq)-control$w):(timePoint-(i*freq)+control$w)])
    }
  }

  # get the parameter for the negative binomial distribution
  # Modification on 13 July 2009 after comment by C. W. Ryan on NAs in the
  # time series
  sumBasevec <- sum(basevec, na.rm=TRUE)
  lengthBasevec <- sum(!is.na(basevec))

  # compute the upper limit of a one sided (1-alpha)*100% prediction interval.
  upPI <- qnbinom(1-control$alpha, sumBasevec + 1/2, (lengthBasevec)/(lengthBasevec + 1))

  # give alarm if the actual value is larger than the upper limit.
  alarm <- observed[timePoint] > upPI

  result <- list(alarm=alarm, upperbound=upPI, disProgObj=disProgObj)
  class(result) = "survRes" # for surveillance system result

  return(result)
}

# 'algo.bayes' calls 'algo.bayesLatestTimepoint' for data points given by range.

algo.bayes <- function(disProgObj, control = list(range = range, b = 0, w = 6, actY = TRUE,alpha=0.05)){

  # Set the default values if not yet set
  if(is.null(control$b)){
    # value from bayes 1
    control$b <- 0
  }
  if(is.null(control$w)){
    # value from bayes 1
    control$w <- 6
  }
  if(is.null(control$alpha)){
    # value from bayes 1
    control$alpha <- 0.05
  }
  if(is.null(control$actY)){
    # value from bayes 1
    control$actY <- TRUE
  }

  # initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)

  count <- 1
  for(i in control$range){
    # call algo.bayesLatestTimepoint
    result <- algo.bayesLatestTimepoint(disProgObj, i, control = control)
    # store the results in the right order
    alarm[count] <- result$alarm
    upperbound[count] <- result$upperbound
    count <- count + 1
  }
  #Add name and data name to control object.
  control$name <- paste("bayes(",control$w,",",control$w*control$actY,",",control$b,")",sep="")
  control$data <- paste(deparse(substitute(disProgObj)))

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj,control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}

algo.bayes1 <- function(disProgObj, control = list(range = range)){
  algo.bayes(disProgObj, control = list(range = control$range, b = 0, w = 6, actY = TRUE,alpha=0.05))
}
algo.bayes2 <- function(disProgObj, control = list(range = range)){
  algo.bayes(disProgObj, control = list(range = control$range, b = 1, w = 6, actY = TRUE,alpha=0.05))
}
algo.bayes3 <- function(disProgObj, control = list(range = range)){
  algo.bayes(disProgObj, control = list(range = control$range, b = 2, w = 4, actY = FALSE,alpha=0.05))
}


