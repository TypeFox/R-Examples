#' Nicely spaced labels along a month
#' 
#' Create dates of certain days of the month for labeling
#' 
#' @return Vector with Dates as returned by \code{\link{as.Date}}.
#' @note Spacing of days is not equal, but set to dertain days of the month! This was originally developed for time series movie frames
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, early 2013
#' @seealso \code{\link{monthAxis}} for automatic determination of npm/npy, \code{\link{as.Date}}, \code{\link{paste}}
#' @keywords chron
#' @export
#' @examples
#' 
#' monthLabs(2014,2014, 3) # 3 days per month
#' monthLabs(2013,2014,  npy=3) # 3 months per year, equally spaced
#' monthLabs(2014,2014,  npy=4) # 4 months per year
#' 
#' # see monthAxis for automatic plot labelling
#' 
#' @param startyear Integer. starting year. DEFAULT: 2002
#' @param stopyear Integer. ending year. DEFAULT: 2018
#' @param npm Integer, one of 1,2,3,6 or 31. Number of labels per month. DEFAULT: 2\cr
#'        npm : days of the month\cr 
#'        1 : first day of each month within the given years\cr 
#'        2 : 1st and 15th day\cr 
#'        3 : 1, 10, 20\cr 
#'        6 : 1, 5, 10, 15, 20, 25.
#'       31 : each day
#' @param npy Integer, one of 1,2,3,4 or 6. Number of labels per year at equally spaced month-beginnings. 
#'        If specified, npm is not considered at all. DEFAULT: NA
#' 
monthLabs <- function(             # returns Date object
startyear=2002,
stopyear=2018,
npm=2,
npy=NA)
{ # start of function
# check for correct input, especially integers:
if(length(startyear)>1) stop("'startyear' has to be one single number.")
if(length(stopyear)>1) stop("'stopyear' has to be one single number.")
if(startyear > stopyear) stop("'startyear' has to be smaller than 'stopyear'.")
if(!abs(startyear - round(startyear)) < .Machine$double.eps^0.5) stop("'startyear' has to be an integer (or close to one).")
if(!abs(stopyear  - round(stopyear )) < .Machine$double.eps^0.5) stop("'stopyear' has to be an integer (or close to one).")
#
# select labels per month or per year
if(is.na(npy)) # the default
{
# check for wrong npm-argument
if(!npm %in% c(1:3,6,31)) stop("wrong 'npm'-value: possible are 1,2,3, 6 or 31.")
# define possible combinations
npmval <- list(npm1=1, npm2=c(1,15), npm3=c(1,10,20), npm4=NA, npm5=NA, npm6=c(1,5,10,15,20,25))
if(npm==31) seq(from=as.Date(paste(startyear,"01-01", sep="-")),
                  to=as.Date(paste(stopyear, "12-31", sep="-")), by="day")
else
# paste the years, months, and (dependent on npm) the days
as.Date(   paste(rep(startyear:stopyear, each=12*npm),
                 rep(1:12,each=npm),
                 npmval[[npm]],       sep="-")
        )
} else # if npy is given:
{
# check for wrong npy-argument
if(!npy %in% c(1:4,6,12)) stop("wrong 'npy'-value: possible are 1,2,3,4,6 or 12.")
# define possible combinations
npyval <- list(npy1=1, npy2=c(1,7), npy3=c(1,5,9), npy4=c(1,4,7,10), npy5=NA, npy6=c(1,3,5,7,9,11))
npyval[[12]] <- 1:12
# paste the years and (dependent on npy) the months
as.Date(paste(rep(startyear:stopyear, each=npy), npyval[[npy]], "01", sep="-"))
} # end of labels per year
} # end of function


# Old stuff in paste command:
                 # numberofyears <- stopyear-startyear ; if(numberofyears<1) numberofyears <- 1
                 # rep(npmval[[npm]], numberofyears*npm),        # old idea, not necessary, as far as I cann tell right now. Maybe paste behaviour changed? Or it was just unnecessary...

