##' @name AnnotateData
##' @export AnnotateData
##' @title Finds midpoints and bands for the within and between cycles.
##' 
##' @description Finds midpoints and bands for the within and between cycles.  This the second of two functions
##' that needs to be called to produce WATS Plots.  \code{AugmentZZZ} is the first.
##' 
##' @param dsLinear The \code{data.frame} to containing the detailed data.
##' @param dvName The name of the dependent/criterion variable.
##' @param centerFunction A function to calculate the center of a subsample.
##' @param spreadFunction A function to calculate the bands of a subsample.
##' @param cycleTallyName The variable name indicating how many cycles have been completed.
##' @param stageIDName The variable name indicating the stage. In a typical interrupted time series, these values are \code{1} before the interruption and \code{2} after.
##' @param stageProgressName The variable name indicating the stage in a decimal form.  This is mostly for internal uses.
##' @param proportionThroughCycleName The variable name indicating how far the point is through a cycle.  For example, 0 degrees would be \code{0}, 180 degrees would be \code{0.5}, 359 degrees would be \code{0.9972}, and 360 degrees would be \code{0}.
##' @param proportionIDName The variable name indicating the ordinal position through a cycle.  
##' @param terminalPointInCycleName The variable name indicating the last point within a given cycle.
##' @return Returns a \code{data.frame} with additional variables <<Say what they are>>.
##' @examples
##' library(Wats)
##' dsLinear <- CountyMonthBirthRate2005Version
##' dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
##' dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
##' 
##' hSpread <- function( scores ) { return( quantile(x=scores, probs=c(.25, .75)) ) }
##' portfolio <- AnnotateData(
##'   dsLinear = dsLinear,
##'   dvName = "BirthRate",
##'   centerFunction = median,
##'   spreadFunction = hSpread
##' )
##' 
##' head(portfolio$dsStageCycle)
##' head(portfolio$dsLinear)
##' head(portfolio$dsPeriodic)

AnnotateData <- function( dsLinear, 
                          dvName,
                          centerFunction,
                          spreadFunction,
                          cycleTallyName="CycleTally", 
                          stageIDName="StageID", 
                          stageProgressName="StageProgress",
                          proportionThroughCycleName="ProportionThroughCycle",
                          proportionIDName="ProportionID",
                          terminalPointInCycleName="TerminalPointInCycle" ) {
  
  pointsInCycle <- max(dsLinear[, proportionIDName])
  testit::assert("The should be at least one point in a cycle", max(pointsInCycle)>=1)
  
  z <- zoo::zooreg(data=dsLinear[, dvName], frequency=pointsInCycle)
  rollingBounds <- zoo::rollapply(data=z, width=pointsInCycle, FUN=spreadFunction)
  
  dsLinear$RollingLower <- NA
  dsLinear$RollingCenter <- NA
  dsLinear$RollingUpper <- NA
  dsLinear$RollingLower[-seq_len(pointsInCycle-1) ] <- rollingBounds[, 1]
  dsLinear$RollingCenter[-seq_len(pointsInCycle-1) ] <- zoo::rollapply(data=z, width=pointsInCycle, FUN=centerFunction)
  dsLinear$RollingUpper[-seq_len(pointsInCycle-1) ] <- rollingBounds[, 2]
  
  
  summarizeStageCycle <- function( d ) {
    positionBounds <- spreadFunction(d[, dvName])
    #   print(positionBounds)
    data.frame(    
      ProportionThroughCycle = mean(d$ProportionThroughCycle, na.rm=TRUE),
      PositionLower = positionBounds[1],
      PositionCenter = centerFunction(d[, dvName]),
      PositionUpper = positionBounds[2]
    )
  }
  dsStageCycle <- plyr::ddply(dsLinear, .variables=c(stageIDName, proportionIDName), .fun=summarizeStageCycle)
  
  dsLinearTemp <- dsLinear[, c("Date", stageIDName, proportionIDName, stageProgressName)]
  colnames(dsLinearTemp)[colnames(dsLinearTemp)==stageIDName] <- "StageIDTime" #Make sure `StageIDTime` matches the two calls below.
  
  dsStageCycleTemp <- dsStageCycle
  colnames(dsStageCycleTemp)[colnames(dsStageCycleTemp)==stageIDName] <- "StageIDBand" #Make sure `StageIDBand` matches the calls below.
  
  dsPeriodic <- merge(x=dsLinearTemp, y=dsStageCycleTemp, by=c(proportionIDName), all.x=TRUE, all.y=TRUE)
  dsPeriodic <- dsPeriodic[order(dsPeriodic[, "Date"], dsPeriodic[, "StageIDTime"], dsPeriodic[, "StageIDBand"]), ]
  
  return( list(dsLinear=dsLinear, dsStageCycle=dsStageCycle, dsPeriodic=dsPeriodic) )
}

# library(Wats)
# dsLinear <- CountyMonthBirthRate2005Version
# dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
# dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
# 
# hSpread <- function( scores ) { return( quantile(x=scores, probs=c(.25, .75)) ) }
# portfolio <- AnnotateData(dsLinear, dvName="BirthRate", centerFunction=median, spreadFunction=hSpread)
# 
# head(portfolio$dsStageCycle)
# head(portfolio$dsLinear)
# head(portfolio$dsPeriodic)
# 
# portfolio <- AnnotateData(dsLinear, dvName="BirthRate", centerFunction=mean, spreadFunction=hSpread)
# 
# head(portfolio$dsStageCycle)
# head(portfolio$dsLinear)
# head(portfolio$dsPeriodic)
