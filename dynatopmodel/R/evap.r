#**********************************************************
# Calculate evapotranspiration curve across part of year
	# dt: time step in hours
	# dStart: start day (1-365)
	# dEnd: end day
	# eMin: evapotranspiration (mm/day) at winter minimum
	# eMax: evapotranspiration (mm/day) at summer maximum
# based on Beven's 1995 Fortran code
#**********************************************************
#load.source("util.r", chdir=T)
require(xts)

###########################################################
# simple potential evaptranspiration calculations
# see Beven (2012) citing Calder et al (1983) who show that this
# is a reasonable approximation to more complex schemes
###########################################################

pEvap <- function(dt=1, dStart=1, numDays=365, eMin=0, eMax=0)
{
  # number days, starting at 0=1st Jan
  days <- (dStart-1):(numDays+dStart-1)

	# ET curve has period of 365, at x=0 and x=364 ET = eMin
	fact <- 1+sin(2*pi*days/365-pi/2)
	# fact goes from 0 at dStart to 2 at day = 365/2
	# calculate daily total ET
	DET <- eMin + 0.5*(eMax-eMin)*fact

	# calculate, start, length and end of all days using sine curve above
  dawn <- 10 - 2.5*fact
  dayLength <- 6 + 4*fact
  sunDown = dawn + dayLength
  dayTotET <- DET

	# if required calculate PE at times steps within day
	# assume no ET outside daylight hours, form is sine curve peaking at
	# midday with total area equal to DET for this day calculated above
  # divide days into stretches of length dt hours

  hrs <- 1:(24%/%dt)*dt

  DET<- NULL
  for(day in 1:numDays)
  {
    # proportion of daylight hours passed
    fract <- (hrs-dawn[day])/dayLength[day]
    # no ET before dawn
    fract[fract<0] <- 0

    # If:
    # T = day total ET
    # l = day length
    # h = hours since dawn
    # f = proportion of day passed since dawn
    # then we have
    # e(h) = 0.5*pi*T/l . sin(pi*h/l) or
    # e(f) = 0.5*pi*T/l . sin(pi*f)
    # total e.t. from dawn to each sample time is
    # 0.5*T*(1-cos(pi*f))
    cE <- -cos(pi*fract)  # omit factor
    # same array shifted right by 1 and cos(0) inserted
    cE2 <- c(-1,cE[1:length(cE)-1])
    # e.t. in interval is difference between cumulative e.t.s
    e <- 0.5*dayTotET[day]*(cE-cE2)

    # deal with values outside daylight hours
    e[which(fract>1)] <- 0
    # add to result
    DET <- c(DET, e)
  }

	return(DET)
}

pe.est <- function(t, pmax=5, t.start=6, t.end=20)
{
	pe <- pmax*sin(pi*(t-t.start)/(t.end-t.start))
	return(pmax(pe, 0))
}

# alternative that takes a start and end local time and returns a time series covering period

# generate an approximated tinme series of evapotranspiration assumming direct relation with sinsuiod insolation
# note daylight times are for mid latitudes
# time step in *hours*
# annual maximum daily total (typically m or mm, be careful)

#' Create sinsuiodal time series of potential evapotranspiration input
#'
#' @details Dynamic TOPMODEL requires a time series of potential
#'   evapotranspiration in order to calculate and remove actual
#'   evapotranpiration from the root zone during a run. Many sophisticated
#'   physical models have been developed for estimating PE and AE, including the
#'   Priestly-Taylor (Priestley and Taylor, 1972) and Penman-Monteith (Montieth,
#'   1965) methods. These, however, require detailed meterological data such as
#'   radiation input and relative humidities that are, in general, difficult to
#'   obtain. Calder (1983) demonstrated that a simple approximation using a
#'   sinusoidal variation in potential evapotranspiration to be a good
#'   approximation to more complex schemes.
#'
#'   If the insolation is also taken to vary sinusoidally through the daylight
#'   hours then, ignoring diurnal meteorological variations, the potential
#'   evapotranspiration during daylight hours for each year day number can be
#'   calculated (for the catchment's latitude). Integration over the daylight
#'   hours allows the daily maximum to be calculated and thus a sub-daily series
#'   generated.
#' @export approx.pe.ts
#' @import xts
#' @importFrom lubridate yday
#' @param dt Time interval in hours
#' @param emin Minimum daily PE total (m or mm)
#' @param emax Maximum daily PE total (m or mm)
#' @param start Start time of returned series in a format that can be coerced into a POSIXct instance. Defaults to start of rainfall data
#' @param end End time for returned series in a format that can be coerced into a POSIXct instance. Defaults to end of rainfall datA
#' @return Time series (xts) of potential evapotranspiration ([L]/[T]) covering the given time range and at the desired interval in m or mm/hr
#' @references Beven, K. J. (2012). Rainfall-runoff modelling : the primer. Chichester, UK, Wiley-Blackwell.
#' @references Calder, I. R. (1986). A stochastic model of rainfall interception. Journal of Hydrology, 89(1), 65-71.
#' @examples
#'\dontrun{
#' # Create PE data for 2012 for use in the Brompton test case
#'
#' require(dynatopmodel)
#'
#' data("brompton")
#'
#'# Generate time series at hourly and 15 minute intervals
#' pe.60 <- approx.pe.ts("2012-01-01", "2012-12-31", dt=1)
#' pe.15 <- approx.pe.ts("2012-01-01", "2012-12-31", dt=0.25)
#'
#' # Check annual totals - should be around 900mm
#' sum(pe.60)*1000
#' sum(pe.15*0.25)*1000
#'
#' # Check maximum daily total on the 1st of July
#' sum(pe.60["2012-07-01"])*1000
#' sum(pe.15["2012-07-01"]*0.25)*1000
#' }
approx.pe.ts <- function(start, end,
                           dt=1,
                           emin=0,
                           emax=5/1000)
{
	# can be supplied as strings
  start <- as.POSIXct(start, origin="1970-01-01")
  end <- as.POSIXct(end, origin="1970-01-01")
  # dt in hours
  # day of year number (uses lubridate) - require whole number of days
  d.start <- lubridate::yday(start)
  d.end <- lubridate::yday(end)

  n.days <- as.double(difftime(end, start, units="days"))
  if(n.days<=0){ stop("Invalid time period specified: start is after end") }

  # use the non xts version
  vals <- pEvap(dt, dStart=d.start, n.days, emin, emax)

  # correct for time step (pe a rate in mm/hr)
  vals <- vals/dt

  # step is in hours = 3600s
  times <- seq(from=start, length.out=length(vals), by=dt*3600)
  len <- min(length(vals), length(times))

  times <- times[1:len]

  res <- xts(vals[1:len], order.by=times)
  return(res)
}
