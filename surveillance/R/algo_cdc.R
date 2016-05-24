###################################################
### chunk number 1: 
###################################################


# Implementation of the CDC surveillance system.
# The system evaluates specified timepoints and gives alarm if it recognizes
# an outbreak for this timepoint.
#

algo.cdcLatestTimepoint <- function(disProgObj, timePoint = NULL, control = list(b = 5, m = 1, alpha=0.025)){

  observed <- disProgObj$observed
  freq <- disProgObj$freq

  # If there is no value in timePoint, then take the last value in observed
  if(is.null(timePoint)){
        timePoint = length(observed)
  }

  # check if the vector observed includes all necessary data.
  if((timePoint-(control$b*freq)-control$m*4) < 1){
        stop("The vector of observed is too short!")
  }

  ######################################################################
  #Find which weeks to take -- hoehle 27.3.2007 - fixed bug taking
  #things in the wrong time order (more recent values)
  ######################################################################
  midx <- seq(-control$m*4-3,control$m*4)
  yidx <- ((-control$b):(-1))*freq 
  baseidx  <- sort(rep(yidx,each=length(midx)) + midx)
  months <- rep(1:((2*control$m+1)*control$b),each=4)
  basevec <- as.integer(by(observed[timePoint + baseidx ],months,sum)) 

  # Create a normal distribution based upper confidence interval
  # (we will use the prediction interval described in 
  # Farrington & Andrew (2003))
  upCi <- mean(basevec)+qnorm(1-control$alpha/2)*sd(basevec)*sqrt(1+1/length(basevec))

  #Counts for the current mounth
  yt0 <- sum(observed[timePoint:(timePoint-3)])
  # Alarm if the actual value is larger than the upper limit.
  alarm <- yt0 > upCi
  # Save aggregated score for later visualisation.
  aggr  <- yt0
  result <- list(alarm=alarm, upperbound=upCi,aggr=aggr)
  class(result) = "survRes" # for surveillance system result
  return(result)
}

# 'algo.cdc' calls 'algo.bayesLatestTimepoint' for data points given by range.

algo.cdc <- function(disProgObj, control = list(range = range, b=5, m=1, alpha=0.025)){
  if(disProgObj$freq != 52) {
    stop("algo.cdc only works for weekly data.")
  }

  # initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  aggr <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)

  #Set control options (standard CDC options)
  if (is.null(control$range)) {
    control$range <- (disProgObj$freq*control$b - control$m):length(disProgObj$observed)
  }
  if (is.null(control$b))        {control$b=5}
  if (is.null(control$m))        {control$m=1} #bug fixed
  if (is.null(control$alpha))    {control$alpha=0.025}

  count <- 1
  for(i in control$range){
    # call algo.cdcLatestTimepoint
    result <- algo.cdcLatestTimepoint(disProgObj, i,control=control)
    # store the results in the right order
    alarm[count] <- result$alarm
    aggr[count] <- result$aggr
    upperbound[count] <- result$upperbound
    count <- count + 1
  }
  #Add name and data name to control object.
  control$name <- paste("cdc(",control$m*4,"*,",0,",",control$b,")",sep="")
  control$data <- paste(deparse(substitute(disProgObj)))

  # Return the vectors-
  # as a special feature CDC objects contain an "aggr" identifier
  # containing the aggregated counts for each week.
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj, control=control, aggr=aggr)

  class(result) = "survRes" # for surveillance system result
  return(result)
}



