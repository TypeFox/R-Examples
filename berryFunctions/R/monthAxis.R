#' Label date axis
#' 
#' Labels date axes at sensible intervals in the time domain of weeks to decades.
#' 
#' @return The dates that were labelled
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2015, update labels and midyear Dec 2015
#' @seealso \code{\link{monthLabs}} for the numbercrunching itself, \code{\link{axis.Date}} with defaults that are less nice.
#' @keywords chron aplot dplot
#' @export
#' @examples
#' 
#' set.seed(007) # for reproducibility
#' Date1 <- sort(as.Date("2013-09-25")+sample(0:150, 30))
#' plot(Date1, cumsum(rnorm(30)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(side=1)
#' monthAxis(1, npm=2, cex.axis=0.5) # fix number of labels per month
#' 
#' plot(Date1, cumsum(rnorm(30)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(labels=FALSE, col.ticks=2)
#' monthAxis(1, format=" ")  # equivalent to axis(labels=FALSE)
#' monthAxis(1)
#' d <- monthAxis(1, labels=letters[1:24], mgp=c(3,2.5,0))
#' d # d covers the full year, thus is longer than n=5
#' 
#' Date2 <- sort(as.Date("2011-07-13")+sample(0:1400, 50))
#' plot(Date2, cumsum(rnorm(50)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(npy=12, format=" ")  # fix number of labels per year
#' monthAxis(tcl=-0.8, lwd.ticks=2, format="%Y/%m", mgp=c(3,1,0))
#' monthAxis(format="", mgp=c(3,2,0)) # International Date format YYYY-mm-dd
#' 
#' plot(Date2, cumsum(rnorm(50)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(midyear=TRUE)
#' abline(v=monthLabs(npm=1), col=8)
#' 
#' Date3 <- sort(as.Date("2011-07-13")+sample(0:1200, 50))
#' plot(Date3, cumsum(rnorm(50)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(1, n=4, font=2)
#' monthAxis(1, col.axis=3) # too many labels with default n=5
#' 
#' # mid-year labels:
#' plot(Date3, cumsum(rnorm(50)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(midyear=TRUE, midargs=list(tcl=-1.2))
#' 
#' # mid-month labels:
#' plot(Date1, cumsum(rnorm(30)), type="l", xaxt="n", ann=FALSE)
#' monthAxis(midmonth=TRUE)
#'
#' # Time axis instead of date axis:
#' plot(as.POSIXct(Sys.time()+c(0,10)*24*3600), 1:2, xaxt="n")
#' monthAxis(n=3)
#' monthAxis()
#'
#' @param side Which \code{\link{axis}} are to be labeled? (can be several). DEFAULT: 1
#' @param timeAxis Logical indicating whether the axis is \code{\link{POSIXct}}, not date. DEFAULT: NA, meaning axis value >1e5
#' @param origin Origin for\code{\link{as.Date}} and \code{\link{as.POSIXct}}. DEFAULT: "1970-01-01"
#' @param startyear Integer. starting year. DEFAULT: NULL = internally computed from \code{\link{par}("usr")}
#' @param stopyear Ditto for ending year. DEFAULT: NULL
#' @param n Approximate number of labels that should be printed (as in \code{\link{pretty}}). DEFAULT: 5
#' @param npm Number of labels per month, overrides n. DEFAULT: NULL = internally computed.
#' @param npy Number of labels per year, overrides npm and n. DEFAULT: NA
#' @param format Format of date, see details in \code{\link{strptime}}. DEFAULT: "\%d.\%m.\\n\%Y"
#' @param labels labels. DEFAULT: format.Date(d, format)
#' @param midyear Place labels in the middle of the year? if TRUE, format default is "\%Y". DEFAULT: FALSE
#' @param midmonth Place labels in the middle of the month? if TRUE, format default is "\%m\\n\%Y". DEFAULT: FALSE
#' @param midargs List of arguments passed to \code{\link{axis}} for the year-start lines without labels. DEFAULT: NULL
#' @param mgp MarGinPlacement, see \code{\link{par}}. The second value is for label distance to axis. DEFAULT: c(3,1.5,0)
#' @param cex.axis CharacterEXpansion (letter size). DEFAULT: 1
#' @param tick Draw tick lines? DEFAULT: TRUE
#' @param las LabelAxisStyle for orientation of labels. DEFAULT: 1 (upright)
#' @param \dots Further arguments passed to \code{\link{axis}}, like \code{lwd, col.ticks, hadj, lty}, ...
#'  
monthAxis <- function(
side=1,
timeAxis=NA,
origin="1970-01-01",
startyear=NULL,
stopyear=NULL,
n=5,
npm=NULL,
npy=NA,
format="%d.%m.\n%Y",
labels=format.Date(d, format),
midyear=FALSE,
midmonth=FALSE,
midargs=NULL,
mgp=c(3,1.5,0),
cex.axis=1,
tick=TRUE,
las=1,
...)
{
# internally needed functions to get Date range from graphic:
getDate <- function(s)
  {
  usr <- par("usr")[if(s%%2) 1:2 else 3:4]
  if(is.na(timeAxis)) timeAxis <- usr[1]>1e5
  if(timeAxis) usr <- as.POSIXct(usr, origin=origin)
  as.Date(usr, origin=origin)
  }
getYear <- function(x) as.numeric(format(x, "%Y"))
# possible combinations of npm, npy:
pos_npm <- c(31,  6,  3,  2,  1, NA, NA,  NA,  NA,  NA)
pos_npy <- c(NA, NA, NA, NA, 12,  6,  4,   3,   2,   1)
pos_dif <- c( 1,  5, 10, 15, 30, 61, 91, 122, 183, 365) # number of days between labels
#
# loop around each side:
for(side_i in side)
  {
  # set from and to:
  dateRange <- getDate(side_i)
  startyear_i <- if(missing(startyear)) getYear(dateRange[1]) else startyear
  stopyear_i  <- if(missing(stopyear )) getYear(dateRange[2]) else stopyear
  # determine npm and npy - desired number of days between labels:
  wish_dif <- as.numeric(diff(dateRange)) / n
  # closest value:
  sel <- which.min(abs(pos_dif - wish_dif))
  # select npm and npy from list
  npm_i <- if(is.null(npm)             ) pos_npm[sel] else npm
  npy_i <- if(is.null(npm) & is.na(npy)) pos_npy[sel] else npy
  # calculate dates
  d <- monthLabs(startyear_i, stopyear_i, npm=npm_i, npy=npy_i)
  # TimeAxis default:
  if(is.na(timeAxis)) timeAxis <- par("usr")[if(side_i%%2) 1 else 3]>1e5
  if(timeAxis) d <- as.POSIXct(d)
  # Label axis
  if(!midyear & ! midmonth) axis(side=side_i, at=d, labels=labels, las=las, mgp=mgp, cex.axis=cex.axis, tick=tick, ...)
  # midyear option:
  if(midyear)
    {
    d <- monthLabs(startyear_i, stopyear_i, npy=2)
    dbor <- d[seq(1,length(d), by=2)] # border dates (=year starting points)
    dmid <- d[seq(2,length(d), by=2)] # mid-year points
    do.call(axis, owa(list(side=side_i, at=dbor, labels=FALSE,las=las, mgp=mgp, 
                           cex.axis=cex.axis, tick=tick), midargs))
    if(missing(format)) format <- "%Y"
    if(missing(mgp)) mgp <- c(3,0.5,0)
    labels <- labels[seq(2,length(d), by=2)]
    axis(side=side_i, at=dmid, labels=labels, las=las, mgp=mgp, cex.axis=cex.axis, tick=FALSE, ...)
    }
  # midmonth option:
  if(midmonth)
    {
    d <- monthLabs(startyear_i, stopyear_i, npm=2)
    dbor <- d[seq(1,length(d), by=2)] # border dates (=year starting points)
    dmid <- d[seq(2,length(d), by=2)] # mid-year points
    do.call(axis, owa(list(side=side_i, at=dbor, labels=FALSE,las=las, mgp=mgp, 
                           cex.axis=cex.axis, tick=tick), midargs))
    if(missing(format)) format <- "%m\n%Y"
    if(missing(mgp)) mgp <- c(3,1.5,0)
    labels <- labels[seq(2,length(d), by=2)]
    axis(side=side_i, at=dmid, labels=labels, las=las, mgp=mgp, cex.axis=cex.axis, tick=FALSE, ...)
    }
  } # End of loop
# output:
return(invisible(d))
} # End of function
