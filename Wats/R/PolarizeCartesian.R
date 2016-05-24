##' @name PolarizeCartesian
##' @export
##' @title Manipulate Cartesian data to use in the WATS polar plot
##' 
##' @description Three operations are performed.  
##' First, within each stage, the first row is repeated at the end, to close the loop.  
##' Second, multiple points are interpolated (still in a Cartesian coordinates) so that the polar graph doesn't have sharp edges.  These sharp edges would be artifacts of the conversion, and not reflect the observed data.
##' Third, the Cartesian points are coverted to polar coordinates.
##' 
##' @param dsLinear The \code{data.frame} to containing the simple linear data.  There should be one record per observation.
##' @param dsStageCycle The \code{data.frame} to containing the reoccurring/periodic bands.  There should be one record per observation per stage.  If there are three stages, this \code{data.frame} should have three times as many rows as \code{dsLinear}.
##' @param yName The variable name containing the dependent/criterion variable.
##' @param stageIDName The variable name indicating which stage the record belongs to.  For example, before the first interruption, the \code{StageID} is \code{1}, and is \code{2} afterwards.
##' @param cycleTallyName The variable name indicating how many \emph{complete} cycles have occurred at that observation.
##' @param proportionThroughCycleName The variable name showing how far through a cycle the observation (or summarized observations) occurred.
##' @param periodicLowerName The variable name showing the lower bound of a stage's periodic estimate.
##' @param periodicCenterName The variable name showing the center estimate of a stage's periodic estimate.
##' @param periodicUpperName The variable name showing the upper bound of a stage's periodic estimate.
##' @param plottedPointCountPerCycle The number of points that are plotted per cycle.  If the polar graph has 'sharp corners', then increase this value.
##' @param graphFloor The value of the criterion/dependent variable at the center of the polar plot.
##' @return Returns a \code{data.frame}.
##' @keywords polar
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
##' rm(dsLinear)
##' 
##' polarized <- PolarizeCartesian(
##'   dsLinear = portfolio$dsLinear, 
##'   dsStageCycle = portfolio$dsStageCycle, 
##'   yName = "BirthRate", 
##'   stageIDName = "StageID"
##' )
##' 
##' library(ggplot2)
##' ggplot(polarized$dsStageCyclePolar, aes(color=factor(StageID))) + 
##'   geom_path(aes(x=PolarLowerX, y=PolarLowerY), linetype=2) +
##'   geom_path(aes(x=PolarCenterX, y=PolarCenterY), size=2) +
##'   geom_path(aes(x=PolarUpperX, y=PolarUpperY), linetype=2) +
##'   geom_path(aes(x=ObservedX, y=ObservedY), data=polarized$dsObservedPolar) +
##'   coord_fixed(ratio=1) +
##'   guides(color=FALSE)

#For a more polished graph, see PolarPeriodic().

PolarizeCartesian <- function(dsLinear, dsStageCycle,
                      yName, stageIDName, 
                      cycleTallyName="CycleTally", 
                      proportionThroughCycleName="ProportionThroughCycle", 
                      periodicLowerName="PositionLower", periodicCenterName="PositionCenter", periodicUpperName="PositionUpper",
                      plottedPointCountPerCycle=120,
                      graphFloor=min(base::pretty(x=dsLinear[, yName]))) {
  #TODO: allow counter-clockwise and arbitrary angle for theta=0
  
  
#   print(dsLinear[, cycleTallyName])
#   print(dsLinear[, proportionThroughCycleName])
#   print(dsLinear[, yName])
  
  closeLoop <- function( d ) {
    d[nrow(d) + 1, ] <- d[1, ] #Within each stage, repeat the first row at the end of the stage's data.frame.
    d[nrow(d), proportionThroughCycleName] <- 1 + d[nrow(d), proportionThroughCycleName]
    return( d )
  }
  interpolateObserved <- function( d, pointsPerCycleCount ) {
    observed <- stats::approx(x = d[, cycleTallyName] + d[, proportionThroughCycleName], 
                              y = d[, yName], 
                              n = pointsPerCycleCount)
    stageProgress <- stats::approx(x = unique(d[, stageIDName]) + 0:1, 
                                   n = pointsPerCycleCount + 1)
    
    base::data.frame(
      ObservedX = observed$x,
      ObservedY = observed$y,
      StageProgress = stageProgress$y[seq_len(pointsPerCycleCount)] #Which chops off the last value.
    )   
  }
  interpolateBand <- function( d, pointsPerCycleCount ) {
    lower <- stats::approx(x=d[, proportionThroughCycleName], y=d[, periodicLowerName], n=pointsPerCycleCount)
    center <- stats::approx(x=d[, proportionThroughCycleName], y=d[, periodicCenterName], n=pointsPerCycleCount)
    upper <- stats::approx(x=d[, proportionThroughCycleName], y=d[, periodicUpperName], n=pointsPerCycleCount)
    
    base::data.frame(
      LowerX = lower$x,
      LowerY = lower$y,
      CenterX = center$x,
      CenterY = center$y,
      UpperX = upper$x,
      UpperY = upper$y
    )   
  }
  polarizeObserved <- function( d, graphFloor=graphFloor ) {
    #After R 3.1.0 has been out for a while, consider using sinpi()`.
    if( nrow(d)==0 ) {
      stageStart <- logical(0)
      stageEnd <- logical(0)
    } else {
      stageStart <- c(TRUE, rep(FALSE, times=nrow(d)-1))
      stageEnd <- c(rep(FALSE, times=nrow(d)-1), TRUE)
    }
    base::data.frame(
      ObservedX = (d$ObservedY - graphFloor) * sin(2 * pi * d$ObservedX),
      ObservedY = (d$ObservedY - graphFloor) * cos(2 * pi * d$ObservedX),
      Theta = pi * 2 * d$ObservedX,
      Radius = d$ObservedY,
      StageProgress = d$StageProgress,
      StageStart = stageStart,
      StageEnd = stageEnd,
      LabelStageStart = ifelse(stageStart, paste0(d$StageID, "S"), ""),
      LabelStageEnd = ifelse(stageEnd, paste0(d$StageID, "E"), ""),
      stringsAsFactors = FALSE      
    )
  }
  polarizeBand <- function( d, graphFloor=graphFloor ) {
    if( nrow(d)==0 ) {
      stageStart <- logical(0)
      stageEnd <- logical(0)
    } else {
      stageStart <- c(TRUE, rep(FALSE, times=nrow(d)-1))
      stageEnd <- c(rep(FALSE, times=nrow(d)-1), TRUE)
    }
    base::data.frame(
      PolarLowerX = (d$LowerY - graphFloor) * sin(2 * pi * d$LowerX),
      PolarLowerY = (d$LowerY - graphFloor) * cos(2 * pi * d$LowerX),  
      PolarCenterX = (d$CenterY - graphFloor) * sin(2 * pi * d$CenterX),
      PolarCenterY = (d$CenterY - graphFloor) * cos(2 * pi * d$CenterX),  
      PolarUpperX = (d$UpperY - graphFloor) * sin(2 * pi * d$UpperX),
      PolarUpperY = (d$UpperY - graphFloor) * cos(2 * pi * d$UpperX),
#       StageProgress = d$StageProgress,
      StageStart = stageStart,
      StageEnd = stageEnd,
      LabelStageStart = ifelse(stageStart, paste0(d$StageID, "S"), ""),
      LabelStageEnd = ifelse(stageEnd, paste0(d$StageID, "E"), ""),
      stringsAsFactors = FALSE
    )
  }
  
  dsObservedInterpolated <- plyr::ddply(dsLinear, .variables=stageIDName, .fun=interpolateObserved, pointsPerCycleCount=plottedPointCountPerCycle)
  dsObservedPolar <- plyr::ddply(dsObservedInterpolated, .variables=stageIDName, .fun=polarizeObserved, graphFloor=graphFloor)

  dsStageCycleClosed <- plyr::ddply(dsStageCycle, .variables=stageIDName, .fun=closeLoop)
  dsStageCycleInterpolated <- plyr::ddply(dsStageCycleClosed, .variables=stageIDName, .fun=interpolateBand, pointsPerCycleCount=plottedPointCountPerCycle)
  dsStageCyclePolar <- plyr::ddply(dsStageCycleInterpolated, .variables=stageIDName, .fun=polarizeBand, graphFloor=graphFloor)
  
  return( list(dsObservedPolar=dsObservedPolar, dsStageCyclePolar=dsStageCyclePolar, GraphFloor=graphFloor) )
}

# library(Wats)
# dsLinear <- CountyMonthBirthRate2005Version
# dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
# dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
# 
# hSpread <- function( scores ) { return( quantile(x=scores, probs=c(.25, .75)) ) }
# portfolio <- AnnotateData(dsLinear, dvName="BirthRate", centerFunction=median, spreadFunction=hSpread)
# rm(dsLinear)
# 
# polarized <- PolarizeCartesian(portfolio$dsLinear, portfolio$dsStageCycle, yName="BirthRate", stageIDName="StageID")
# 
# library(ggplot2)
# ggplot(polarized$dsStageCyclePolar, aes(color=factor(StageID))) + 
#   geom_path(aes(x=PolarLowerX, y=PolarLowerY), linetype=2) +
#   geom_path(aes(x=PolarCenterX, y=PolarCenterY), size=2) +
#   geom_path(aes(x=PolarUpperX, y=PolarUpperY), linetype=2) +
#   geom_path(aes(x=ObservedX, y=ObservedY), data=polarized$dsObservedPolar) +
#   coord_fixed(ratio=1) +
#   guides(color=FALSE)
# 
# #For a more polished graph, see PolarPeriodic().
