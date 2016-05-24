##############################################################################################
# Calculates life table estimator, Survival function and standard deviation without covariates
##############################################################################################

# Description
# Constructs life table and estimates hazard rate, survival function and corresponding standard errors

# Input
# dataSet: Original data set of format data.frame in short format
# timeColumn: Scalar character giving the column name of the time variable
# censColumn: Scalar character giving the column name of the censoring variable

# Output
# data frame with columns
  # Rows: Each row of the table is an half open interval with width=1 (character)
  # n: Number of individuals at risk in a given time interval (integer)
  # events: Observed number of events in a given time interval (integer)
  # dropouts: Observed number of dropouts in a given time interval (integer)
  # atRisk: Estimated number of individuals at risk, corrected by dropouts
  # hazard: Estimated risk of death (without covariates) in a given time interval
  # seHazard: Estimated standard deviation of estimated hazard
  # S: Estimated survival curve
  # seS: Estimated standard deviation of estimated survival function

lifeTable <- function (dataSet, timeColumn, censColumn, intervalBorders=NULL) {
  
  # Input Checks
  if(!is.data.frame(dataSet)) {stop("Data set is not of type data.frame! Please convert the data!")}
  if(!(length(timeColumn)==1 & is.character(timeColumn))) {stop("The column name is not correctly specified! The argument must be a scalar as character.")}
  if(!(length(censColumn)==1 & is.character(censColumn))) {stop("The column name is not correctly specified! The argument must be a scalar as character.")}
  if(!(is.null (intervalBorders) )) {
    if(!(is.character(intervalBorders))) {stop("The interval borders are not in character format! Please give the appropriate format, e. g. [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q) with a_q beeing the number of intervals")}
    MaxTime <- max(dataSet [, timeColumn])
    if(!(length(intervalBorders)==MaxTime)) {stop("The interval borders have not the same length as the number of unique intervals! Please give the appropriate format, e. g. [0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q) with a_q beeing the number of intervals")}
  }
  
  # Data transformation to life table
  atRiskInitial <- dim(dataSet) [1]
  dataSet <- dataSet [order(dataSet [, timeColumn]), ]
  formulaInput <- as.formula(paste(censColumn, timeColumn, sep="~"))
  events <- aggregate(formula=formulaInput, data=dataSet, FUN=function (x) sum(x)) [, 2]
  dropouts <- aggregate(formula=formulaInput, data=dataSet, FUN=function (x) sum(1-x)) [, 2]
  atRiskInput <- c(atRiskInitial, atRiskInitial-cumsum(events+dropouts))
  atRiskInput <- atRiskInput [-length(atRiskInput)]
  times <- as.numeric(names(table(as.numeric(as.character(dataSet [, timeColumn])))))
  Index <- which(diff(times)>1)
  while(any(diff(times)>1)) {
    Index <- which(diff(times)>1) [1]
    atRiskInput <- c(atRiskInput [1:Index], 
                     atRiskInput [Index] - (events[Index] + dropouts[Index]), 
                     atRiskInput [(Index+1):length(atRiskInput)])
    events <- c(events [1:Index], 0, events [(Index+1):length(events)])
    dropouts <- c(dropouts [1:Index], 0, dropouts [(Index+1):length(dropouts)])
    times <- c(times [1:Index], times[Index] + 1, times [(Index+1):length(times)])
  }
  
  # Correction of first observed category:
  # It is assumed that no event and no dropouts occur in unobserved intervals
  # Each line must correspond to one interval
  if(times[1]!=1) {
    atRiskInput <- c(rep(atRiskInput [1], times[1] -1), atRiskInput)
    events <- c(rep(0, times[1] - 1), events)
    dropouts <- c(rep(0, times[1] - 1), dropouts)
    times <- c(1:(times[1] - 1), times)
  }
  
  # Estimation of hazards, survival function cumulative hazard and standard deviations
  atRisk <- atRiskInput - dropouts/2
  haz <- events / atRisk
  S <- cumprod(1 - haz)
  sehaz <- sqrt((haz - haz^2) / atRisk)
  seS <- S * sqrt(cumsum(haz / (1 - haz) / (atRisk)))
  cumHazard <- cumsum(haz)
  seCumHazard <- sqrt(cumsum(events / atRisk^2))
  
  # Construct interval borders
  if(is.null(intervalBorders)) {
    RowNames <- paste("[", c(0, times [-length(times)]), ", ", times, ")", sep = "")
  }
  else {
    RowNames <- intervalBorders
  }

  # RowNames <- RowNames [-length(RowNames)]
  Output <- data.frame(n=atRiskInput, events = events, dropouts = dropouts, atRisk = atRisk,
             hazard = haz, seHazard = sehaz,
             S = S, seS = seS, 
             cumHazard=cumHazard, seCumHazard=seCumHazard, 
             row.names = RowNames)
  
  # Exclude last row because theoretically the hazard will be 1 because all subjects die. 
  # However the estimator can be zero if no event occurs. Therefore it is not reliable
  # Output <- Output [-dim(Output) [1], ]
  Output <- list(Output=Output)
  class(Output) <- "discSurvLifeTable"
  return(Output)
}

print.discSurvLifeTable <- function (x, ...) {
  x <- x [[1]]
  for(i in 1:dim(x) [2]) {
    x [, i] <- round(x [, i], 4)
  }
  
  print(x)
}