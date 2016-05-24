##' @name PolarPeriodic
##' @export
##' @title Polar Plot with Periodic Elements
##' 
##' @description Shows the interrupted time series in Cartesian coordinates and its a periodic/cyclic components.
##' 
##' @param dsLinear The \code{data.frame} to containing the simple linear data.  There should be one record per observation.
##' @param dsStageCyclePolar The \code{data.frame} to containing the bands for a single period.  There should be one record per theta per stage.  If there are three stages, this \code{data.frame} should have three times as many rows as \code{dsLinear}.
##' @param xName The variable name containing the date.
##' @param yName The variable name containing the dependent/criterion variable.
##' @param stageIDName The variable name indicating which stage the record belongs to.  For example, before the first interruption, the \code{StageID} is \code{1}, and is \code{2} afterwards.
# 
##' @param periodicLowerName The variable name showing the lower bound of a stage's periodic estimate.
##' @param periodicUpperName The variable name showing the upper bound of a stage's periodic estimate.
##' @param paletteDark A vector of colors used for the dark/heavy graphical elements.  The vector should have one color for each \code{StageID} value.  If no vector is specified, a default will be chosen, based on the number of stages.
##' @param paletteLight A vector of colors used for the light graphical elements.  The vector should have one color for each \code{StageID} value.  If no vector is specified, a default will be chosen, based on the number of stages.
##' @param changePoints A vector of values indicate the interruptions between stages.  It typically works best as a \code{Date} or a \code{POSIXct} class.
##' @param changePointLabels The text plotted above each interruption.
##' @param drawObservedLine A boolean value indicating if the longitudinal observed line should be plotted (whose values are take from \code{dsLinear}).
##' @param drawPeriodicBand A boolean value indicating if the bands should be plotted (whose values are take from the \code{periodicLowerName} and \code{periodicUpperName} fields).
##' @param drawStageLabels A boolean value indicating if the stage labels should be plotted (whose values are take from \code{dsLinear}).
##' @param drawRadiusLabels A boolean value indicating if the gridline/radius labels should be plotted (whose values are take from \code{tickLocations}).
##' @param jaggedPointSize The size of the observed data points.
##' @param jaggedLineSize The size of the line connecting the observed data points.
##' 
##' @param bandAlphaDark The amount of transparency of the band appropriate for a stage's \emph{x} values.
##' @param bandAlphaLight The amount of transparency of the band comparison stages for a given \emph{x} value.
##' @param changeLineAlpha The amount of transparency marking each interruption.
##' @param colorLabels The color for \code{cardinalLabels} and \code{originLabel}.
##' @param colorGridlines The color for the gridlines.
##' @param labelColor The color of the text labels imposed on the line.
##' @param changeLineSize The width of a line marking an interruption.
##' @param tickLocations The desired locations for ticks showing the value of the criterion/dependent variable.
##' @param graphFloor The value of the criterion/dependent variable at the center of the polar plot.
##' @param graphCeiling The value of the criterion/dependent variable at the outside of the polar plot.
##' 
##' @param cardinalLabels The four labels placed  where `North', `East', `South', and `West' typically are.
##' @param originLabel Explains what the criterion variable's value is at the origin.  Use \code{NULL} if no explanation is desired.
##' @param plotMargins A vector of four \code{numeric} values, specifying the number of lines in the bottom, left, top and right margins.
##' 
##' @return Returns a \code{grid} graphical object (ie, a \href{http://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.grob.html}{\code{grob}}.)
##' @keywords polar
##' @examples
##' library(grid)
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
##'   portfolio$dsLinear, 
##'   portfolio$dsStageCycle, 
##'   yName = "BirthRate", 
##'   stageIDName = "StageID"
##' )
##' 
##' grid.newpage()
##' PolarPeriodic(
##'   dsLinear = polarized$dsObservedPolar, 
##'   dsStageCyclePolar = polarized$dsStageCyclePolar, 
##'   yName = "Radius", 
##'   stageIDName = "StageID", 
##'   cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1")
##' )
##' 
##' grid.newpage()
##' PolarPeriodic(
##'   dsLinear = polarized$dsObservedPolar, 
##'   dsStageCyclePolar = polarized$dsStageCyclePolar, 
##'   yName = "Radius", 
##'   stageIDName = "StageID", 
##'   drawPeriodicBand = FALSE
##' )
##' 
##' grid.newpage()
##' PolarPeriodic(
##'   dsLinear = polarized$dsObservedPolar, 
##'   dsStageCyclePolar = polarized$dsStageCyclePolar, 
##'   yName = "Radius", 
##'   stageIDName = "StageID", 
##'   drawObservedLine = FALSE, 
##'   cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1")
##' )

PolarPeriodic <- function(dsLinear, dsStageCyclePolar,
                          xName, yName, stageIDName, 
                          periodicLowerName = "PositionLower", periodicUpperName = "PositionUpper",
                          paletteDark = NULL, paletteLight = NULL, 
                          changePoints = NULL, changePointLabels = NULL,
                          drawObservedLine = TRUE, drawPeriodicBand = TRUE, 
                          drawStageLabels = FALSE, drawRadiusLabels = FALSE,
                          jaggedPointSize = 2, jaggedLineSize = 1, 
                          bandAlphaDark = .4, bandAlphaLight = .15, 
                          colorLabels = "gray50", colorGridlines = "gray80", labelColor="orange3",
                          changeLineAlpha = .5, changeLineSize = 3,
                          tickLocations = base::pretty(x = dsLinear[, yName]),
                          graphFloor = min(tickLocations),
                          graphCeiling = max(tickLocations),
                          cardinalLabels = NULL, originLabel = paste0("The origin represents ", graphFloor, ";\nthe perimeter represents ", graphCeiling, "."), 
                          plotMargins = c(3.5, 2, .5, 2)
                          ) {
  
  testit::assert("The `dsStageCyclePolar` must have a valid column called `PolarLowerX`.  Typically this is generated by `Wats::PolarizeCartesian()`.", !is.null(dsStageCyclePolar$PolarLowerX))
  testit::assert("The `dsStageCyclePolar` must have a valid column called `PolarLowerY`.  Typically this is generated by `Wats::PolarizeCartesian()`.", !is.null(dsStageCyclePolar$PolarLowerY))
  testit::assert("The `dsStageCyclePolar` must have a valid column called `PolarUpperX`.  Typically this is generated by `Wats::PolarizeCartesian()`.", !is.null(dsStageCyclePolar$PolarUpperX))
  testit::assert("The `dsStageCyclePolar` must have a valid column called `PolarUpperY`.  Typically this is generated by `Wats::PolarizeCartesian()`.", !is.null(dsStageCyclePolar$PolarUpperY))
  
  tickLocationsPolar <- tickLocations - min(tickLocations)
  
  graphRadius <- graphCeiling - graphFloor
  vpRange <- c(-graphRadius, graphRadius) * 1.02
  stages <- base::sort(base::unique(dsLinear[, stageIDName]))
  stageCount <- length(stages)
  #     testit::assert("The number of unique `StageID` values should be 1 greater than the number of `changePoints`.", stageCount==1+length(changePoints))
  if( !is.null(changePoints) ) testit::assert("The number of `changePoints` should equal the number of `changeLabels`.", length(changePoints)==length(changePointLabels))
  if( !is.null(paletteDark) ) testit::assert("The number of `paletteDark` colors should equal the number of unique `StageID` values.", stageCount==length(paletteDark))
  if( !is.null(paletteLight) ) testit::assert("The number of `paletteLight` colors should equal the number of unique `StageID` values.", stageCount==length(paletteLight))

  if( is.null(paletteDark) ) {
    if( length(stages) <= 4L) paletteDark <- RColorBrewer::brewer.pal(n=10L, name="Paired")[c(2L,4L,6L,8L)] #There's not a risk of defining more colors than levels
    else paletteDark <- colorspace::rainbow_hcl(n=length(stages), l=40)
  }  
  if( is.null(paletteLight) ) {
    if( length(stages) <= 4L) paletteLight <- RColorBrewer::brewer.pal(n=10L, name="Paired")[c(1L,3L,5L,7L)] #There's not a risk of defining more colors than levels
    else paletteLight <- colorspace::rainbow_hcl(n=length(stages), l=70)
  }   
#   grid.rect() #For exploring nested viewports
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow=1, ncol=1, respect=T), gp=grid::gpar(cex=0.6, fill=NA)))
#   grid.rect() #For exploring nested viewports
  grid::pushViewport(grid::viewport(layout.pos.col=1, layout.pos.row=1)) #This simple viewport is very important for the respected aspect ratio of 1.
#   grid.rect() #For exploring nested viewports
  grid::grid.text(originLabel, x=0, y=0, hjust=-.1, vjust=-.2, gp=grid::gpar(cex=1.5, col=colorLabels, lineheight=.8), default.units="npc")
  grid::pushViewport(grid::plotViewport(margins=plotMargins)) 
#   grid.rect() #For exploring nested viewports
  grid::pushViewport(grid::dataViewport(xscale=vpRange, yscale=vpRange, name="plotRegion"))
#   grid.rect() #For exploring nested viewports

  grid::grid.lines(x=c(-graphRadius,graphRadius), y=c(0,0), gp=grid::gpar(col=colorGridlines, lty=3), default.units="native")
  grid::grid.lines(x=c(0,0), y=c(-graphRadius,graphRadius), gp=grid::gpar(col=colorGridlines, lty=3), default.units="native")
  grid::grid.circle(x=0, y=0, r=tickLocationsPolar, default.units="native", gp=grid::gpar(col=colorGridlines))
  if( drawRadiusLabels ) {
#     grid::grid.text(tickLocations, x=tickLocationsPolar, y=0, default.units="native", 
#                     gp=grid::gpar(col=colorGridlines), just=c(-.1, 1.1))
    grid::grid.text(tickLocations, x=tickLocationsPolar/sqrt(2), y=-tickLocationsPolar/sqrt(2), default.units="native", 
                    gp=grid::gpar(col=colorLabels), just=c(-.05, 1.05))
  }
  grid::grid.text(cardinalLabels, x=c(0, graphRadius, 0, -graphRadius), y=c(graphRadius, 0, -graphRadius, 0), gp=grid::gpar(cex=2, col=colorLabels), default.units="native")
    
#   lg <- grid::polylineGrob(x=dsStageCyclePolar$PolarLowerX, y=dsStageCyclePolar$PolarLowerY, id=dsStageCyclePolar$StageID, gp=grid::gpar(col=paletteDark, lwd=2), default.units="native", name="l") #summary(lg) #lg$gp
#   grid::grid.draw(lg)   
#   cg <- grid::polylineGrob(x=dsStageCyclePolar$PolarCenterX, y=dsStageCyclePolar$PolarCenterY, id=dsStageCyclePolar$StageID, gp=grid::gpar(col=paletteDark, lwd=2), default.units="native", name="l") #summary(lg) #lg$gp
#   grid::grid.draw(cg)    
#   ug <- grid::polylineGrob(x=dsStageCyclePolar$PolarUpperX, y=dsStageCyclePolar$PolarUpperY, id=dsStageCyclePolar$StageID, gp=grid::gpar(col=paletteDark, lwd=2), default.units="native", name="l") #summary(lg) #lg$gp
#   grid::grid.draw(ug)

  if( drawPeriodicBand ) {
    for( stageID in stages ) {
      lowerX <- dsStageCyclePolar[dsStageCyclePolar$StageID==stageID, "PolarLowerX"]
      lowerY <- dsStageCyclePolar[dsStageCyclePolar$StageID==stageID, "PolarLowerY"]
      upperX <- dsStageCyclePolar[dsStageCyclePolar$StageID==stageID, "PolarUpperX"]
      upperY <- dsStageCyclePolar[dsStageCyclePolar$StageID==stageID, "PolarUpperY"]  
      
      x <- c(lowerX, rev(upperX))
      y <- c(lowerY, rev(upperY))
      grid::grid.polygon(x=x, y=y, default.units="native", gp=grid::gpar(fill=paletteDark[stageID], col="transparent", alpha=bandAlphaDark))
    }
  }
  
  if( drawObservedLine ) {
#     gObserved <- grid::polylineGrob(x=dsLinear$ObservedX, y=dsLinear$ObservedY, id=dsLinear$StageID, 
#                                     gp=grid::gpar(col=paletteDark, lwd=jaggedLineSize), 
#                                     default.units="native", name="l") 
#     grid::grid.draw(gObserved)
    for( stage in stages) {
      dsStage <- dsLinear[stage<=dsLinear$StageProgress & dsLinear$StageProgress<=(stage+1), ]
      
      gObserved <- grid::polylineGrob(x=dsStage$ObservedX, y=dsStage$ObservedY, 
                                      gp=grid::gpar(col=paletteDark[stage], lwd=jaggedLineSize), 
                                      default.units="native", name="l") 
      grid::grid.draw(gObserved)
      
    }
  }

  if( drawStageLabels ) {
    gLabelStart <- grid::textGrob(label=dsLinear$LabelStageStart, x=dsLinear$ObservedX, y=dsLinear$ObservedY, 
                                  gp=grid::gpar(col=labelColor, lwd=jaggedLineSize), 
                                  default.units="native", name="l") 
    grid::grid.draw(gLabelStart)
    gLabelEnd <- grid::textGrob(label=dsLinear$LabelStageEnd, x=dsLinear$ObservedX, y=dsLinear$ObservedY, 
                                gp=grid::gpar(col=labelColor, lwd=jaggedLineSize), 
                                default.units="native", name="l") 
    grid::grid.draw(gLabelEnd)
  }
  grid::upViewport(n=4)
}

# library(grid)
# library(Wats)
# dsLinear <- CountyMonthBirthRate2005Version
# dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
# dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
# 
# hSpread <- function( scores ) { return( quantile(x=scores, probs=c(.25, .75)) ) }
# portfolio <- AnnotateData(dsLinear, dvName="BirthRate", centerFunction=median, spreadFunction=hSpread)
# rm(dsLinear)
# 
# polarized <- PolarizeCartesian(portfolio$dsLinear, portfolio$dsStageCycle, yName="BirthRate", stageIDName="StageID", plottedPointCountPerCycle=3600)
# 
# grid.newpage()
# PolarPeriodic(dsLinear=polarized$dsObservedPolar, polarized$dsStageCyclePolar, drawRadiusLabels=TRUE, drawStageLabels=TRUE, yName="Radius", stageIDName="StageID", drawPeriodicBand=FALSE)

# grid.newpage()
# PolarPeriodic(dsLinear=polarized$dsObservedPolar, polarized$dsStageCyclePolar, yName="Radius", stageIDName="StageID", drawPeriodicBand=FALSE)

# grid.newpage()
# PolarPeriodic(dsLinear=polarized$dsObservedPolar, polarized$dsStageCyclePolar, yName="Radius", stageIDName="StageID", cardinalLabels=c("Jan1", "Apr1", "July1", "Oct1"))

# # 
# grid.newpage()
# PolarPeriodic(dsLinear=polarized$dsObservedPolar, polarized$dsStageCyclePolar, yName="Radius", stageIDName="StageID", drawObservedLine=FALSE, cardinalLabels=c("Jan1", "Apr1", "July1", "Oct1"))
