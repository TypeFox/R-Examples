## ----set_root_directory, echo=FALSE, results='hide'----------------------
#It works better if the root directory is set in its own chunk.
library(knitr)
knitr::opts_knit$set(root.dir = "../")

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
    out.width="600px", #This affects only the markdown, not the png.  Height will scale appropriately.
    fig.path='figure_okc_fertility_intercensal_rmd/',
    dpi=200
)

base::options(width=120) #So the output is 50% wider than the default.
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

## ----CartesianRolling, fig.height=4.8---------------------------------------------------------------------------------
# dsLinear <- utils::read.csv("./Datasets/CountyMonthBirthRate2014Version.csv", stringsAsFactors=FALSE)
# dsLinear$Date <- base::as.Date(dsLinear$Date)
# dsLinear <- dsLinear[dsLinear$CountyName=="oklahoma", ] 

# Uncomment this line to use the version built into the package.  By default, it uses the
# CSV to promote reproducible research, since the CSV format is more open and accessible to more software.
dsLinearAll <- CountyMonthBirthRate2014Version
dsLinear <- dsLinearAll[dsLinearAll$CountyName=="oklahoma", ] 

dsLinear <- Wats::AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")

portfolioCartesian <- Wats::AnnotateData(dsLinear, dvName="BirthRate", centerFunction=stats::median, spreadFunction=hSpread)

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

middlePanel <- CartesianRolling(
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

## ----CartesianPeriodic------------------------------------------------------------------------------------------------
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
cartesianPeriodic <- cartesianPeriodic + xScale + darkTheme 
print(cartesianPeriodic)

## ----PolarPeriodic, fig.height=6.5*2/3--------------------------------------------------------------------------------
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
    widths = grid::unit(c(1,1), c("null", "null")), 
    heights = grid::unit(c(1,.5), c("null", "null"))
  ), 
  gp = grid::gpar(cex=1, fill=NA)
))

## Create top left panel
grid::pushViewport(grid::viewport(layout.pos.col=1, layout.pos.row=1))
topLeftPanel <- Wats::PolarPeriodic(  
  dsLinear = portfolioPolar$dsObservedPolar, 
  dsStageCyclePolar = portfolioPolar$dsStageCyclePolar, 
  yName = "Radius", 
  stageIDName = "StageID",
  cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1")
)
grid::upViewport()

## Create top right panel
grid::pushViewport(grid::viewport(layout.pos.col=2, layout.pos.row=1))
topRighttPanel <- Wats::PolarPeriodic(
  dsLinear = portfolioPolar$dsObservedPolar, 
  dsStageCyclePolar = portfolioPolar$dsStageCyclePolar, 
  yName = "Radius", 
  stageIDName = "StageID",
  drawObservedLine = FALSE,
  cardinalLabels = c("Jan1", "Apr1", "July1", "Oct1"), 
  originLabel = NULL
)
grid::upViewport()

## Create bottom panel
grid::pushViewport(grid::viewport(layout.pos.col=1:2, layout.pos.row=2, gp=grid::gpar(cex=1)))
print(cartesianPeriodic, vp=vpLayout(x=1:2, y=2)) #Print across both columns of the bottom row.
upViewport()

## ----ConfirmatoryFrequentist, fig.height=6.5--------------------------------------------------------------------------
dsLinear <- Wats::AugmentYearDataWithMonthResolution(dsLinear=dsLinear, dateName="Date")

tsData <- stats::ts(
  data = dsLinear$BirthRate, 
  start = as.integer(dsLinear[1, c("Year", "Month")]), 
  end = as.integer(dsLinear[nrow(dsLinear), c("Year", "Month")]),
  frequency = 12L
)

#Create unsmoothed and smoothed version
seasonalClassic <- stats::decompose(tsData)
plot(seasonalClassic)

#Watch out, the 2nd & 3rd columns have swapped positions, compared to `decompose()`
seasonalLoess <- stats::stl(x = tsData, s.window = "periodic") 
plot(seasonalLoess)

# Seasonality is accounted for without a smoother
lag1 <- 1L #Significant for many different values of lag, including 1
y <- seasonalClassic$trend[(lag1+1):length(seasonalClassic$trend)]
y1 <- seasonalClassic$trend[1:(length(seasonalClassic$trend)-lag1)]
# step <- c(rep(0L, times=sum(dsLinear$StageID==1L)-lag1), rep(1L, times=sum(dsLinear$StageID==2L)))
step <- dsLinear$StageID[(lag1+1):length(seasonalClassic$trend)] - 1L
dsClassic <- data.frame(y=y, y1=y1, step=step)
rm(lag1, y, y1, step)
fitClassic <-  glm(y ~ 1 + step + y1, data=dsClassic)
summary(fitClassic)

#Seasonality is accounted for after a loess is fit through it.
lag1 <- 1L #Significant for many different values of lag, including 1
trendLineLoess <- as.numeric(seasonalLoess$time.series[, 2])
y <- trendLineLoess[(lag1+1):length(trendLineLoess)]
y1 <- trendLineLoess[1:(length(trendLineLoess) - lag1)]
# step <- c(rep(0L, times=sum(dsLinear$StageID==1L)-lag1), rep(1L, times=sum(dsLinear$StageID==2L)))
step <- dsLinear$StageID[(lag1+1):length(trendLineLoess)] - 1L
dsLoess <- data.frame(y=y, y1=y1, step=step)
rm(lag1, y, y1, step)
fitLoess <-  glm(y ~ 1 + step + y1, data=dsLoess)
summary(fitLoess)

## ----ConfirmatoryBayesFactors, fig.height=6.5-------------------------------------------------------------------------
# Seasonality is accounted for without a smoother
beforeClassic <- seasonalClassic$trend[dsLinear$StageID==1L]
afterClassic <- seasonalClassic$trend[dsLinear$StageID==2L]

#Set the number of MCMC iterations
mcmcRepCount <- 1000#000
#Determine if the intermediate progress should be displayed
showMcmcProgress <- FALSE

(gClassic <- BayesSingleSub::trendtest.Gibbs.AR(beforeClassic[!is.na(beforeClassic)], afterClassic[!is.na(afterClassic)], return.chains=F, iterations=mcmcRepCount, progress=showMcmcProgress))
(mcClassic <- BayesSingleSub::trendtest.MC.AR(beforeClassic[!is.na(beforeClassic)], afterClassic[!is.na(afterClassic)], iterations=mcmcRepCount, progress=showMcmcProgress))
# summary(mcClassic)
# coda::gelman.diag(g$chains) #it needs multiple chains to asses.

#Seasonality is accounted for after a loess is fit through it.
beforeLoess <- seasonalLoess$time.series[dsLinear$StageID==1L, 2]
afterLoess <- seasonalLoess$time.series[dsLinear$StageID==2L, 2]
BayesSingleSub::trendtest.Gibbs.AR(beforeLoess, afterLoess, iterations=mcmcRepCount, progress=showMcmcProgress)
BayesSingleSub::trendtest.MC.AR(beforeLoess, afterLoess, iterations=mcmcRepCount, progress=showMcmcProgress)


## ----session_info, echo=FALSE-----------------------------------------------------------------------------------------
base::cat("Report created by", base::Sys.info()["user"], "at", base::strftime(base::Sys.time(), "%c, %z"))
utils::sessionInfo()

