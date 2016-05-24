######################################################################
# Takes a data frame with dates of individual
# cases and create an aggregated sts time series object for these
# data with aggregation occuring at the desired scale.
#
# Parameters:
#  linelist - a data frame containing individual case information, one per line
#  dateCol - a character string denoting the column name in case containing
#            the relevant date variable to aggregate
#  aggregate.by - aggregation block length given as a string compatible with
#       seq.Date -- see \link{seq.Date} for further details.
#
# Author: Michael Hoehle
# Date LaMo: 04 Jan 2014
######################################################################

linelist2sts <- function(linelist,dateCol,aggregate.by=c("1 day", "1 week", "7 day", "1 week", "1 month", "3 month", "1 year"),dRange=NULL,
                      epochInPeriodStr=switch(aggregate.by, "1 day"="1","1 week"="%u", "1 month"="%d","3 month"="%q","1 year"="%j"),
                      startYearFormat=switch(aggregate.by,"1 day"="%Y","7 day"="%G","1 week"="%G","1 month"="%Y","3 month"="%Y","1 year"="%Y"),
                      startEpochFormat=switch(aggregate.by,"1 day"="%j","7 day"="%V","1 week"="%V","1 month"="%m","3 month"="%Q","1 year"="1")
                      ) {

  ##Check aggregate.by argument
  aggregate.by <- match.arg(aggregate.by, c("1 day", "1 week", "7 day", "1 week", "1 month", "3 month", "1 year"))
  
  #If no dRange let it be the range of the dateCol
  if (is.null(dRange)) {
    dRange <- range(linelist[,dateCol],na.rm=TRUE)
  }
  if (aggregate.by != "1 day") {
    ##Move dates back to first of each epoch unit
    dRange <-  dRange - as.numeric(formatDate(dRange,epochInPeriodStr)) + 1
  }
  
  #Add exactly one time step to dRange to ensure that cut
  #contains the last level as well. We use 'seq' to ensure
  #that even weeks/days with no data are present in the factor.
  maxDate <- seq(max(dRange),length.out=2,by=aggregate.by)[-1]
  dates <- seq(min(dRange), maxDate, by=aggregate.by)

  #Make a table containing the specific number of cases. Note that this
  #needs to occur using a cut statement
  lvl <- cut(linelist[,dateCol], breaks=dates,right=FALSE)

  observed <- table(lvl)
  epoch <- as.Date(names(observed))

  #Translate "by" to freq string
  freq <- switch(aggregate.by,"1 day"=365,"7 day"=52,"1 week"=52,"1 month"=12,"3 month"=4,"1 year"=1)

  startYear <- as.numeric(formatDate(min(dates),startYearFormat))
  startEpoch <- as.numeric(formatDate(min(dates),startEpochFormat))
                  
  observed <- matrix(observed,ncol=1)

  #Create S4 object
  sts <- new("sts",epoch=as.numeric(epoch),observed=observed, alarm=0*observed, epochAsDate=TRUE,freq=freq,start=c(startYear,startEpoch))

  #Return
  return(sts)
}
