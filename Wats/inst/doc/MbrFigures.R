## ----set_root_directory, echo=FALSE, results='hide'----------------------
#It works better if the root directory is set in its own chunk.
# library(knitr)
# knitr::opts_knit$set(root.dir = "../")

## ----set_options, echo=FALSE, results='hide'--------------------------------------------------------------------------
library(base)
library(utils)
library(stats)
library(grDevices)
library(grid)
library(plyr) 
library(scales)
library(ggplot2) 
library(boot) 
library(Wats)

knitr::opts_chunk$set(
    comment=NA, 
    tidy=FALSE,
    fig.width=6.5, 
    fig.height=1.6, 
    out.width="600px", #This affects only the markdown, not the underlying png file.  The height will be scaled appropriately.
    fig.path='figure_mbr_rmd/',
    # dev = "pdf", #Uncomment to produce vector graphics for publication; however png (the default) works better within html (ie, the format of this vignette).
    dpi=400
)

base::options(width=120) #So the output is 50% wider than the default.
pdf.options(useDingbats = FALSE) #Otherwise, the circles don't get plotted correctly in some graphs
if( base::Sys.info()["sysname"] == "Windows" )
  grDevices::windows.options(antialias = "cleartype")

## ----Definitions------------------------------------------------------------------------------------------------------
changeMonth <- base::as.Date("1996-02-15") #as.Date("1995-04-19") + lubridate::weeks(39) = "1996-01-17"

vpLayout <- function(x, y) { grid::viewport(layout.pos.row=x, layout.pos.col=y) }

fullSpread <- function( scores ) { 
  return( base::range(scores) ) #A new function isn't necessary.  It's defined in order to be consistent.
}
hSpread <- function( scores ) { 
  return( stats::quantile(x=scores, probs=c(.25, .75)) ) 
}
seSpread <- function( scores ) { 
  return( base::mean(scores) + base::c(-1, 1) * stats::sd(scores) / base::sqrt(base::sum(!base::is.na(scores))) ) 
}
bootSpread <- function( scores, conf=.68 ) {
  plugin <- function( d, i ) { base::mean(d[i]) }

  distribution <- boot::boot(data=scores, plugin, R=99) #999 for the publication
  ci <- boot::boot.ci(distribution, type = c("bca"), conf=conf)
  return( ci$bca[4:5] ) #The fourth & fifth elements correspond to the lower & upper bound.
}

darkTheme <- ggplot2::theme(
  axis.title          = ggplot2::element_text(color="gray30", size=9),
  axis.text.x         = ggplot2::element_text(color="gray30", hjust=0),
  axis.text.y         = ggplot2::element_text(color="gray30"),
  axis.ticks.length   = grid::unit(0, "cm"),
  axis.ticks.margin   = grid::unit(.00001, "cm"),
#   panel.grid.minor.y  = element_line(color="gray95", size=.1),
#   panel.grid.major    = element_line(color="gray90", size=.1),
  panel.margin        = grid::unit(c(0, 0, 0, 0), "cm"),
  plot.margin         = grid::unit(c(0, 0, 0, 0), "cm")
)
lightTheme <- darkTheme + ggplot2::theme(
  axis.title          = ggplot2::element_text(color="gray80", size=9),
  axis.text.x         = ggplot2::element_text(color="gray80", hjust=0),
  axis.text.y         = ggplot2::element_text(color="gray80"),
  panel.grid.minor.y  = ggplot2::element_line(color="gray99", size=.1),
  panel.grid.major    = ggplot2::element_line(color="gray95", size=.1)
)
dateSequence <- base::seq.Date(from=base::as.Date("1990-01-01"), to=base::as.Date("1999-01-01"), by="years")
xScale       <- ggplot2::scale_x_date(breaks=dateSequence, labels=scales::date_format("%Y"))
xScaleBlank  <- ggplot2::scale_x_date(breaks=dateSequence, labels=NULL) #This keeps things proportional down the three frames.

## ----Figure2IndividualBasic-------------------------------------------------------------------------------------------
# dsLinearAll <- utils::read.csv("./Datasets/CountyMonthBirthRate2005Version.csv", stringsAsFactors=FALSE)
# dsLinearAll$Date <- base::as.Date(dsLinearAll$Date)
# dsLinearOkc <- dsLinearAll[dsLinearAll$CountyName=="oklahoma", ] 

# Uncomment the next two lines to use the version built into the package.  By default, it uses the
# CSV to promote reproducible research, since the CSV format is more open and accessible to more software.
dsLinearAll <- CountyMonthBirthRate2005Version
dsLinearOkc <- dsLinearAll[dsLinearAll$CountyName=="oklahoma", ]

dsLinearOkc <- Wats::AugmentYearDataWithMonthResolution(dsLinear=dsLinearOkc, dateName="Date")

portfolioCartesian <- Wats::AnnotateData(dsLinearOkc, dvName="BirthRate", centerFunction=stats::median, spreadFunction=hSpread)

Wats::CartesianRolling(
  dsLinear = portfolioCartesian$dsLinear, 
  xName = "Date", 
  yName = "BirthRate", 
  stageIDName = "StageID", 
  changePoints = changeMonth, 
  changePointLabels = "Bombing Effect"
)

## ----Figure2Stylized, fig.height=4.8----------------------------------------------------------------------------------
topPanel <- Wats::CartesianRolling(
  dsLinear = portfolioCartesian$dsLinear, 
  xName = "Date", 
  yName = "BirthRate", 
  stageIDName = "StageID", 
  changePoints = changeMonth, 
  yTitle = "General Fertility Rate",
  changePointLabels = "Bombing Effect", 
  drawRollingBand = FALSE, 
  drawRollingLine = FALSE
)

middlePanel <- Wats::CartesianRolling(
  dsLinear = portfolioCartesian$dsLinear, 
  xName = "Date", 
  yName = "BirthRate", 
  stageIDName = "StageID", 
  changePoints = changeMonth, 
  yTitle = "General Fertility Rate",
  changePointLabels = "", 
  drawRollingBand = FALSE, 
  drawJaggedLine = FALSE
)

bottomPanel <- Wats::CartesianRolling(
  dsLinear = portfolioCartesian$dsLinear, 
  xName = "Date", 
  yName = "BirthRate", 
  stageIDName = "StageID", 
  changePoints = changeMonth, 
  yTitle = "General Fertility Rate", 
  changePointLabels = "", 
#   drawRollingBand = FALSE, 
  drawJaggedLine = FALSE
)

topPanel <- topPanel + xScale + darkTheme
middlePanel <- middlePanel + xScale + darkTheme
bottomPanel <- bottomPanel + xScaleBlank + darkTheme

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1)))
print(topPanel, vp=vpLayout(1, 1))
print(middlePanel, vp=vpLayout(2, 1))
print(bottomPanel, vp=vpLayout(3, 1))
grid::popViewport()

## ----Figure4Basic-----------------------------------------------------------------------------------------------------
cartesianPeriodic <- Wats::CartesianPeriodic(
  portfolioCartesian$dsLinear, 
  portfolioCartesian$dsPeriodic, 
  xName = "Date", 
  yName = "BirthRate",
  stageIDName = "StageID", 
  changePoints = changeMonth, 
  changePointLabels = "Bombing Effect",
  yTitle = "General Fertility Rate",
  drawPeriodicBand = TRUE #The only difference from the simple linear graph above
)
print(cartesianPeriodic)

## ----Figure4Stylized--------------------------------------------------------------------------------------------------
cartesianPeriodic <- cartesianPeriodic + xScale + darkTheme 
print(cartesianPeriodic)

## ----Figure5, fig.height=3, fig.width=3, out.width="300px"------------------------------------------------------------
portfolioPolar <- Wats::PolarizeCartesian( 
  dsLinear = portfolioCartesian$dsLinear, 
  dsStageCycle = portfolioCartesian$dsStageCycle, 
  yName = "BirthRate", 
  stageIDName = "StageID", 
  plottedPointCountPerCycle = 7200
)

grid::grid.newpage()
Wats::PolarPeriodic(
  dsLinear = portfolioPolar$dsObservedPolar, 
  dsStageCycle = portfolioPolar$dsStageCyclePolar, 
  yName = "Radius",
  stageIDName = "StageID", 
  drawPeriodicBand = FALSE,
  drawStageLabels = TRUE, 
  drawRadiusLabels = TRUE, 
  cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1")
)


## ----Figure6, fig.height=6.5*2/3--------------------------------------------------------------------------------------
portfolioPolar <- Wats::PolarizeCartesian(
  dsLinear = portfolioCartesian$dsLinear, 
  dsStageCycle = portfolioCartesian$dsStageCycle, 
  yName = "BirthRate", 
  stageIDName = "StageID", 
  plottedPointCountPerCycle = 7200
)

grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout=grid::grid.layout(
    nrow = 2, ncol = 2, respect = TRUE, 
    widths = unit(c(1, 1), c("null", "null")), 
    heights = unit(c(1, .5), c("null", "null"))
  ), 
  gp = grid::gpar(cex=1, fill=NA)
))

## Create top left panel
grid::pushViewport(grid::viewport(layout.pos.col=1, layout.pos.row=1))
topLeftPanel <- Wats::PolarPeriodic(  
  dsLinear = portfolioPolar$dsObservedPolar, 
  dsStageCyclePolar = portfolioPolar$dsStageCyclePolar, 
  yName = "Radius", 
  stageIDName = "StageID", #graphCeiling=7, 
  cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1")
)
grid::upViewport()

## Create top right panel
grid::pushViewport(grid::viewport(layout.pos.col=2, layout.pos.row=1))
topRighttPanel <- Wats::PolarPeriodic(
  dsLinear = portfolioPolar$dsObservedPolar, 
  dsStageCyclePolar = portfolioPolar$dsStageCyclePolar, 
  yName = "Radius", 
  stageIDName = "StageID", #graphCeiling=7, 
  drawObservedLine = FALSE,
  cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1"), 
  originLabel = NULL
)
grid::upViewport()

## Create bottom panel
grid::pushViewport(grid::viewport(layout.pos.col=1:2, layout.pos.row=2, gp=grid::gpar(cex=1)))
print(cartesianPeriodic, vp=vpLayout(x=1:2, y=2)) #Print across both columns of the bottom row.
grid::upViewport()

## ----Figure7, fig.height=6.5------------------------------------------------------------------------------------------
# dsLinearAll <- Wats::AugmentYearDataWithMonthResolution(dsLinear=CountyMonthBirthRate2005Version, dateName="Date")

#Identify the average size of the fecund population
plyr::ddply(dsLinearAll, "CountyName", plyr::summarize, Mean=base::mean(FecundPopulation))


GraphRowComparison <- function( rowLabel="", countyName="oklahoma", spreadFunction=hSpread, changeMonth=as.Date("1996-02-15") ) {
  dsLinear <- dsLinearAll[dsLinearAll$CountyName==countyName, ]
  dsLinear <- Wats::AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")
  portfolioCartesian <- Wats::AnnotateData(dsLinear, dvName="BirthRate", centerFunction=stats::median, spreadFunction=spreadFunction)
  portfolioPolar <- Wats::PolarizeCartesian(dsLinear=portfolioCartesian$dsLinear, dsStageCycle=portfolioCartesian$dsStageCycle, yName="BirthRate", stageIDName="StageID", plottedPointCountPerCycle=7200)
  cartesianPeriodic <- Wats::CartesianPeriodic(portfolioCartesian$dsLinear, portfolioCartesian$dsPeriodic, xName="Date", yName="BirthRate", stageIDName="StageID", changePoints=changeMonth, changePointLabels=""  )
  
  grid::pushViewport(grid::viewport(
    layout=grid::grid.layout(nrow=1, ncol=3, respect=F, widths=grid::unit(c(1.5,1,3), c("line", "null", "null"))), 
    gp=grid::gpar(cex=1, fill=NA)
  ))
  grid::pushViewport(grid::viewport(layout.pos.col=1))
  grid.rect(gp=gpar(fill="gray90", col=NA))
  grid.text(rowLabel, rot=90)
  grid::popViewport()
  
  grid::pushViewport(grid::viewport(layout.pos.col=2))
  polarPeriodic <- Wats::PolarPeriodic(dsLinear=portfolioPolar$dsObservedPolar, dsStageCyclePolar=portfolioPolar$dsStageCyclePolar, drawObservedLine=FALSE, yName="Radius", stageIDName="StageID", originLabel=NULL, plotMargins=c(0,0,0,0))
  grid::popViewport()
  
  grid::pushViewport(grid::viewport(layout.pos.col=3))
  print(cartesianPeriodic + xScale + lightTheme, vp=vpLayout(x=1, y=1))
  grid::popViewport()
  grid::popViewport() #Finish the row
}

counties <- c("comanche", "cleveland", "oklahoma", "tulsa", "rogers")
countyNames <- c("Comanche", "Cleveland", "Oklahoma", "Tulsa", "Rogers")

grid.newpage()
grid::pushViewport(grid::viewport(layout=grid.layout(nrow=length(counties), ncol=1), gp=grid::gpar(cex=1, fill=NA)))
for( i in base::seq_along(counties) ) {
  grid::pushViewport(grid::viewport(layout.pos.row=i))
  GraphRowComparison(countyName=counties[i], rowLabel=countyNames[i])
  grid::popViewport()
}
grid::popViewport()

## ----Figure7AllCounties, fig.height=6.5 * 12/5------------------------------------------------------------------------
counties <- base::sort(base::unique(dsLinearAll$CountyName))
countyNames <- c("Canadian", "Cleveland", "Comanche", "Creek", "Logan", "McClain", "Oklahoma", "Osage", "Pottawatomie", "Rogers", "Tulsa", "Wagoner")

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout=grid.layout(nrow=base::length(counties), ncol=1), gp=grid::gpar(cex=1, fill=NA)))
for( i in base::seq_along(counties) ) {
  grid::pushViewport(grid::viewport(layout.pos.row=i))
  GraphRowComparison(countyName=counties[i], rowLabel=countyNames[i])
  grid::popViewport()
}
grid::popViewport()

## ----Figure8, fig.height=6.5 * 4/5------------------------------------------------------------------------------------
spreads <- c("hSpread", "fullSpread", "seSpread", "bootSpread")
spreadNames <- c("H-Spread", "Range", "+/-1 SE", "Bootstrap")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow=base::length(spreads), ncol=1), gp=grid::gpar(cex=1, fill=NA)))
for( i in base::seq_along(spreads) ) {
  grid::pushViewport(grid::viewport(layout.pos.row=i))
  GraphRowComparison(spreadFunction=base::get(spreads[i]), rowLabel=spreadNames[i])
  grid::upViewport()
}
grid::upViewport()

## ----session_info, echo=FALSE-----------------------------------------------------------------------------------------
base::cat("Report created by", base::Sys.info()["user"], "at", base::strftime(base::Sys.time(), "%c, %z"))
utils::sessionInfo()

