### R code from vignette source 'Rnw/algo_rki.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: algo_rki.Rnw:96-214
###################################################


# Implementation of the Robert-Koch Institute (RKI) surveillance system.
# The system evaluates specified timepoints and gives alarm if it recognizes
# an outbreak for this timepoint.
#
# Features:
# Choice between the different RKI sub-systems (difference in reference values).

algo.rkiLatestTimepoint <- function(disProgObj, timePoint = NULL, control = list(b = 2, w = 4, actY = FALSE)){

  observed <- disProgObj$observed
  freq <- disProgObj$freq

  # If there is no value in timePoint, then take the last value in observed
  if(is.null(timePoint)){
        timePoint = length(observed)
  }

  # check if the vector observed includes all necessary data.
  if((timePoint-(control$b*freq)-control$w) < 1){
        stop("The vector of observed is too short!")
  }

  # Extract the reference values from the historic time series
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

  # compute the mean.
  mu <- mean(basevec)

  if(mu > 20){ # use the normal distribution.
    # comupte the standard deviation.
    sigma <- sqrt(var(basevec))
    # compute the upper limit of the 95% CI.
    upCi <- mu + 2 * sigma
  }
  else{ # use the poisson distribution.
    # take the upper limit of the 95% CI from the table CIdata.txt.
    #data("CIdata", envir=environment())   # only local assignment -> SM: however, should not use data() here
    #CIdata <- read.table(system.file("data", "CIdata.txt", package="surveillance"), header=TRUE)
    #SM: still better: use R/sysdata.rda (internal datasets being lazy-loaded into the namespace environment)
    # for the table-lookup mu must be rounded down.
    mu <- floor(mu)
    # we need the third column in the row mu + 1
    upCi <- CIdata[mu + 1, 3]
  }
  # give alarm if the actual value is larger than the upper limit.
  alarm <- observed[timePoint] > upCi

  result <- list(alarm=alarm, upperbound=upCi)
  class(result) = "survRes" # for surveillance system result
  return(result)
}

# 'algo.rki' calls 'algo.bayesLatestTimepoint' for data points given by range.

algo.rki <- function(disProgObj, control = list(range = range, b = 2, w = 4, actY = FALSE)){
  # Set the default values if not yet set
  if(is.null(control$b)){
    # value from rki 3
    control$b <- 2
  }
  if(is.null(control$w)){
    # value from rki 3
    control$w <- 4
  }
  if(is.null(control$actY)){
    # value from rki 3
    control$actY <- FALSE
  }

  # initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)

  count <- 1
  for(i in control$range){
    #hoehle Debug:
    #print(i)
    # call algo.rki1LatestTimepoint
    result <- algo.rkiLatestTimepoint(disProgObj, i, control = control)
    # store the results in the right order
    alarm[count] <- result$alarm
    upperbound[count] <- result$upperbound
    count <- count + 1
  }

  #Add name and data name to control object.
  control$name <- paste("rki(",control$w,",",control$w*control$actY,",",control$b,")",sep="")
  control$data <- paste(deparse(substitute(disProgObj)))

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj, control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}

algo.rki1 <- function(disProgObj, control = list(range = range)) {
  algo.rki(disProgObj, control = list(range = control$range, b = 0, w = 6, actY = TRUE))
}
algo.rki2 <- function(disProgObj, control = list(range = range)){
  algo.rki(disProgObj, control = list(range = control$range, b = 1, w = 6, actY = TRUE))
}
algo.rki3 <- function(disProgObj, control = list(range = range)){
  algo.rki(disProgObj, control = list(range = control$range, b = 2, w = 4, actY = FALSE))
}



