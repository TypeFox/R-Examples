##' @name CartesianPeriodic
##' @export
##' @title Linear Plot with Periodic Elements
##' 
##' @description Shows the interrupted time series in Cartesian coordinates and its a periodic/cyclic components.
##' 
##' @param dsLinear The \code{data.frame} to containing the simple linear data.  There should be one record per observation.
##' @param dsPeriodic The \code{data.frame} to containing the reoccurring/periodic bands.  There should be one record per observation per stage.  If there are three stages, this \code{data.frame} should have three times as many rows as \code{dsLinear}.
##' @param xName The variable name containing the date.
##' @param yName The variable name containing the dependent/criterion variable.
##' @param stageIDName The variable name indicating which stage the record belongs to.  For example, before the first interruption, the \code{StageID} is \code{1}, and is \code{2} afterwards.
##' @param periodicLowerName The variable name showing the lower bound of a stage's periodic estimate.
##' @param periodicUpperName The variable name showing the upper bound of a stage's periodic estimate.
##' @param paletteDark A vector of colors used for the dark/heavy graphical elements.  The vector should have one color for each \code{StageID} value.  If no vector is specified, a default will be chosen, based on the number of stages.
##' @param paletteLight A vector of colors used for the light graphical elements.  The vector should have one color for each \code{StageID} value.  If no vector is specified, a default will be chosen, based on the number of stages.
##' @param changePoints A vector of values indicate the interruptions between stages.  It typically works best as a \code{Date} or a \code{POSIXct} class.
##' @param changePointLabels The text plotted above each interruption.
##' @param drawPeriodicBand A boolean value indicating if the bands should be plotted (whose values are take from the \code{periodicLowerName} and \code{periodicUpperName}.

##' @param jaggedPointSize The size of the observed data points.
##' @param jaggedLineSize The size of the line connecting the observed data points.
##' 
##' @param bandAlphaDark The amount of transparency of the band appropriate for a stage's \emph{x} values.
##' @param bandAlphaLight The amount of transparency of the band comparison stages for a given \emph{x} value.
##' @param changeLineAlpha The amount of transparency marking each interruption.
##' @param changeLineSize The width of a line marking an interruption.
##' 
##' @param title The string describing the plot.
##' @param xTitle The string describing the \emph{x}-axis.
##' @param yTitle The string describing the \emph{y}-axis. 
##' 
##' @return Returns a \code{ggplot2} graphing object
##' @keywords Cartesian
##' @examples
##' library(Wats) #Load the package
##' changeMonth <- base::as.Date("1996-02-15")
##' dsLinear <- CountyMonthBirthRate2005Version
##' dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ]
##' dsLinear <- AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
##' hSpread <- function( scores ) { return( quantile(x=scores, probs=c(.25, .75)) ) }
##' portfolio <- AnnotateData(
##'     dsLinear,
##'     dvName = "BirthRate",
##'     centerFunction = median,
##'     spreadFunction = hSpread
##' )
##' 
##' CartesianPeriodic(
##'   portfolio$dsLinear,
##'   portfolio$dsPeriodic,
##'   xName = "Date",
##'   yName = "BirthRate",
##'   stageIDName = "StageID",
##'   changePoints = changeMonth,
##'   changePointLabels = "Bombing Effect"
##' )



CartesianPeriodic <- function(dsLinear, dsPeriodic,
                              xName, yName, stageIDName, 
                              periodicLowerName="PositionLower", periodicUpperName="PositionUpper",
                              paletteDark=NULL, paletteLight=NULL, 
                              changePoints=NULL, changePointLabels=NULL,
                              drawPeriodicBand=TRUE,
                              jaggedPointSize=2, jaggedLineSize=.5, 
                              bandAlphaDark=.4, bandAlphaLight=.15, 
                              changeLineAlpha=.5, changeLineSize=3,
                              title=NULL, xTitle=NULL, yTitle=NULL ) {
  
  stages <- base::sort(base::unique(dsLinear[, stageIDName]))
  stageCount <- length(stages)
  testit::assert("The number of unique `StageID` values should be 1 greater than the number of `changePoints`.", stageCount==1+length(changePoints))
  if( !is.null(changePoints) ) testit::assert("The number of `changePoints` should equal the number of `changeLabels`.", length(changePoints)==length(changePointLabels))
  if( !is.null(paletteDark) ) testit::assert("The number of `paletteDark` colors should equal the number of unique `StageID` values.", stageCount==length(paletteDark))
  if( !is.null(paletteLight) ) testit::assert("The number of `paletteLight` colors should equal the number of unique `StageID` values.", stageCount==length(paletteLight))
  
  p <- ggplot2::ggplot(dsLinear, ggplot2::aes_string(x=xName, y=yName))
  
  if( is.null(paletteDark) ) {
    if( length(stages) <= 4L) paletteDark <- RColorBrewer::brewer.pal(n=10, name="Paired")[c(2,4,6,8)] #There's not a risk of defining more colors than levels
    else paletteDark <- colorspace::rainbow_hcl(n=length(stages), l=40)
  }  
  if( is.null(paletteLight) ) {
    if( length(stages) <= 4L) paletteLight <- RColorBrewer::brewer.pal(n=10, name="Paired")[c(1,3,5,7)] #There's not a risk of defining more colors than levels
    else paletteLight <- colorspace::rainbow_hcl(n=length(stages), l=70)
  }  
    
  for( stage in stages ) {
    dsStageLinear <- dsLinear[stage<=dsLinear$StageProgress & dsLinear$StageProgress<=(stage+1), ]
    
    if( drawPeriodicBand ) {
      for( stageInner in stages ) {
        dsStagePeriodic <- dsPeriodic[(stage<=dsPeriodic$StageProgress) & (dsPeriodic$StageProgress<=(stage+1)) & (dsPeriodic$StageIDBand==stageInner), ]
        ribbonAlpha <- ifelse(stage==stageInner, bandAlphaDark, bandAlphaLight)
        #p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin=periodicLowerName, ymax=periodicUpperName, y=NULL), data=dsStagePeriodic, 
        #                     fill=paletteDark[stageInner], color=NA, alpha=ribbonAlpha, na.rm=T)
        
        p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(y=periodicLowerName, ymin=periodicLowerName, ymax=periodicUpperName), data=dsStagePeriodic, 
                              fill=paletteDark[stageInner], color=NA, alpha=ribbonAlpha, na.rm=T)
      }
    }
        
    p <- p + ggplot2::geom_line(size=jaggedLineSize, color=paletteDark[stage], data=dsStageLinear)    
    p <- p + ggplot2::geom_point(shape=1, color=paletteLight[stage], data=dsStageLinear, size=jaggedPointSize)
  }

  if( !is.null(changePoints) ) {
    for( i in seq_along(changePoints) )  {
      p <- p + ggplot2::geom_vline(xintercept=as.integer(changePoints[i]), color=paletteLight[i+1], alpha=changeLineAlpha, size=changeLineSize)
      p <- p + ggplot2::annotate("text", x=changePoints[i], y=Inf, vjust=1.1, color=paletteLight[i+1], label=changePointLabels[i])
    }
  }

  p <- p + ggplot2::theme_minimal()
  p <- p + ggplot2::theme(legend.position="none")
  p <- p + ggplot2::labs(title=title, x=xTitle, y=yTitle)
  
  return( p )
}

# dsLinear <- CountyMonthBirthRate2005Version
# dsLinear[dsLinear$CountyName=="oklahoma", ]
# dsLinear <- Wats::AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
# 
# hSpread <- function( scores ) { return( quantile(x=scores, probs=c(.25, .75)) ) }
# portfolio <- Wats::AnnotateData(dsLinear, dvName="BirthRate", centerFunction=median, spreadFunction=hSpread)
# 
# CartesianPeriodic(portfolio$dsLinear, portfolio$dsPeriodic, xName="Date", yName="BirthRate", stageIDName="StageID", changePoints=changeMonth, changePointLabels="Bombing Effect",
#                    drawPeriodicBand=FALSE)
# CartesianPeriodic(portfolio$dsLinear, portfolio$dsPeriodic, xName="Date", yName="BirthRate", stageIDName="StageID", changePoints=changeMonth, changePointLabels="Bombing Effect")
