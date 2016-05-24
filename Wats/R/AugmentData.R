##' @name AugmentCycleData
##' @aliases AugmentYearDataWithMonthResolution AugmentYearDataWithSecondResolution
##' @export AugmentYearDataWithMonthResolution AugmentYearDataWithSecondResolution
##' @usage AugmentYearDataWithMonthResolution( dsLinear, dateName ) 
##' AugmentYearDataWithSecondResolution( dsLinear, dateName ) 
##' 
##' @title Calculates variables necessary for WATS Plots
##' 
##' @description Calculates variables necessary for WATS Plots.  This the first of two functions
##' that needs to be called to produce WATS Plots.  \code{AnnotateData} is the second.
##' 
##' @param dsLinear The \code{data.frame} to containing the detailed data.
##' @param dateName The variable name in \code{dsLinear} containing the date or datetime value.
## @param stageIDName The variable name indicating the stage. In a typical interrupted time series, these values are \code{1} before the interruption and \code{2} after.
##' @return Returns a \code{data.frame} with additional variables: \code{CycleTally}, \code{ProportionThroughCycle}, \code{ProportionID}, and \code{TerminalPointInCycle}.
##' @examples
##' library(Wats)
##' dsLinear <- CountyMonthBirthRate2005Version
##' dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
##' dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
##' head(dsLinear)
##' 
AugmentYearDataWithMonthResolution <- function( dsLinear, dateName ) {
  yearOfEvent <- lubridate::year(dsLinear[, dateName])

  minYearOfEvent <- base::min(yearOfEvent)
  dsLinear$CycleTally <- (yearOfEvent - minYearOfEvent)
  monthsThroughTheYear <- lubridate::month(dsLinear[, dateName]) - .5
  monthsInTheYear <- 12L
  dsLinear$ProportionThroughCycle <- monthsThroughTheYear /  monthsInTheYear
  dsLinear$ProportionID <- base::rank(dsLinear$ProportionThroughCycle, ties.method="max") / base::max(dsLinear$CycleTally + 1)
  dsLinear$StartingPointInCycle <- (dsLinear$ProportionID==base::min(dsLinear$ProportionID))
  dsLinear$TerminalPointInCycle <- (dsLinear$ProportionID==base::max(dsLinear$ProportionID))
  
  SummarizeWithinStage <- function( d ) {
    isMin <- (base::min(d[, dateName]) < d[, dateName])
    return( d$StageID + isMin*0.5 )
  }
  dsLinear$StageProgress <- base::unlist(plyr::dlply(dsLinear, "StageID", SummarizeWithinStage))
  return( dsLinear )
}
AugmentYearDataWithSecondResolution <- function( dsLinear, dateName ) {
  yearOfEvent <- lubridate::year(dsLinear[, dateName])
  firstOfYear <- base::ISOdate(year=yearOfEvent, month=1, day=1, tz="GMT")
  lastOfYear <- firstOfYear + lubridate::years(1)  #ISOdate(year=yearOfEvent + 1, month=1, day=1, tz="GMT") 
  
  minYearOfEvent <- min(yearOfEvent)
  dsLinear$CycleTally <- (yearOfEvent - minYearOfEvent)
  secondsThroughTheYear <- base::as.integer(base::difftime(time1=dsLinear[, dateName], firstOfYear, units="sec")) - .5
  secondsInTheYear <- base::as.integer(base::difftime(lastOfYear, firstOfYear, units="sec"))
  dsLinear$ProportionThroughCycle <- secondsThroughTheYear /  secondsInTheYear
  
  SummarizeWithinCycle <- function( d ) {
    d$ProportionID <- base::rank(d$ProportionThroughCycle, ties.method="max")
    d$StartingPointInCycle <- (d$ProportionID==base::min(d$ProportionID))
    d$TerminalPointInCycle <- (d$ProportionID==base::max(d$ProportionID)) 
    return( d )
  }
  dsLinear <- plyr::ddply(dsLinear, .variables="CycleTally", SummarizeWithinCycle) #base::transform,
#                           ProportionID)
  
  #dsLinear$ProportionID <- as.integer(round(rank(dsLinear$ProportionThroughCycle, ties.method="max") / max(dsLinear$CycleTally + 1)))
#   dsLinear$ProportionID <- rank(dsLinear$ProportionThroughCycle, ties.method="max") / max(dsLinear$CycleTally + 1)
#   dsLinear$StartingPointInCycle <- (dsLinear$ProportionID==min(dsLinear$ProportionID))
#   dsLinear$TerminalPointInCycle <- (dsLinear$ProportionID==max(dsLinear$ProportionID))  
#   dsLinear <- plyr::ddply(dsLinear, 
#                     "CycleTally", 
#                     transform, 
#                     TerminalPointInCycle=(rank(ProportionThroughCycle)==max(rank(ProportionThroughCycle))))  
  SummarizeWithinStage <- function( d ) {
    #     minValue <- min(d[, dateName])
    #     maxValue <- max(d[, dateName])
    #     isBetween <- ( (min(d[, dateName]) < d[, dateName]) & (d[, dateName] < max(d[, dateName])))
    isMin <-  (base::min(d[, dateName]) < d[, dateName])
    return( d$StageID + isMin*0.5 )
  }
  dsLinear$StageProgress <- base::unlist(plyr::dlply(dsLinear, "StageID", SummarizeWithinStage))
#   dsLinear$StageProgress <- plyr::daply(dsLinear, "StageID", SummarizeWithinStage)
  return( dsLinear )
}

# library(Wats)
# dsLinear <- CountyMonthBirthRate2005Version
# dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
# # dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
# dsLinear
# 
# dsLinear$Date <- as.POSIXct(dsLinear$Date, tz="GMT")
# dsLinear <- AugmentYearDataWithSecondResolution(dsLinear=dsLinear, dateName="Date")
#   
