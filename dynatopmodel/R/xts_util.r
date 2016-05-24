#------------------------------------------------------------------------------------------------------------
# Time series utilities

require(zoo)
require(xts)

#' Resample observation data at a new time interval
#'
#' @description Takes a list of time series and resample to a new interval.
#' @details Time series of observation data are often of different temporal resolutions, however the input to most hydrological models, as is the case with the Dynamic TOPMODEL, requires those data at the same interval. This provides a method to resample a collection of such data to a single interval.
#' @export aggregate_obs
#' @param obs List of times series (zoo) objects with a POSIXct index.
#' @param dt New time interval, hours.
#' @param is.rate If TRUE then these are rates i.e m/hr. Otherwise they are absolute values across the interval and are scaled before return by a factor equal to the ratio of the old interval to the new interval.
#' @return The list of observations resampled at the new interval.
#' @examples
#' # Resample Brompton rainfall and PE data to 15 minute intervals
#' require(dynatopmodel)
#' data("brompton")
#'
#' obs <- aggregate_obs(list("rain"=brompton$rain, "pe"=brompton$pe), dt=15/60)
#'
#' # check totals for Sept - Oct 2012
#' sum(obs$rain*15/60, na.rm=TRUE)
#' sum(brompton$rain, na.rm=TRUE)
#'

aggregate_obs <- function(obs, dt, is.rate=TRUE)
{
  #nms <- names(obs)
  obs <- sapply(obs, disaggregate_xts, dt=dt, is.rate=is.rate)
#  names(obs)<-nms
  return(obs)
}

# useful POSIX const (sys dependent?)
t.origin <- function(){as.POSIXct("1970-01-01 00:00:00")}

# get the entries in the xts / zoo series at the given time period
subset_zoo <- function(obj, start, end)
{
  start <- as.POSIXct(start)
  end <- as.POSIXct(end)
  # 	require(intervals)
  #   int1 <-
  #   int2 <- Intervals(c(start, end))
  if(!is.zoo(obj))
  {
    warning("Non-zoo object supplied to subset_zoo")
    return(obj)
  }
  set1 <- which(index(obj) >= start)
  set2 <- which(index(obj) <= end)
  #  if(length(set1)==0){return(obj[])}
  ind <- intersect(set1,set2)
  return(obj[ind])
}

# initialise an xts object from a vectors of date times
def.xts <- function(tms, nser=1, val=0)
{
  dat <- matrix(val, ncol=nser, nrow=length(tms))
  xts(dat, order.by=tms)
}
subset.xts <- function(ser, start=NULL, end=NULL)
{

  if(is.null(ser)){return(ser)}
  if(!is(ser, "xts")){stop("xts object required")}
  if(is.null(start)){start <- start(ser)}
  if(is.null(start)){start <- end(ser)}

  start <- as.POSIXct(start)
  end <- as.POSIXct(end)
  sel <- paste0 (start, "::", end)
  return(ser[sel])
}


# split supplied series into another with higher resolution dt or supply factor
# optionally supply time range to subset input series
disaggregate_xts <- function(ser,
                             ser.start=start(ser),
                             ser.end=end(ser),
                             fact=NULL,    # factor by which to divide existing time interval e.g fact=4 converts houly intervals to quarter hours
                             dt=NULL,
                             is.rate=TRUE)  # IF true then each value in teh series is replicated across the finer time intervals, otherwise the initial value is averaged averaged over
{
  # if the set is NULL then return
  if(is.null(ser)){return(ser)}
  if(!is.zoo(ser)){stop("Time series required")}
  if(is.null(fact)& is.null(dt)){stop("Supply either aggregation factor or new time interval")}
  ser <- subset_zoo(ser, ser.start, ser.end)
  tms <- as.double(index(ser))
  dt.ser <- round(mean(diff(tms)/3600,na.rm=T),2)

  # frequency of the input series to the
  if(is.null(fact))
  {
    # disaggregation factor (nearest)
    tryCatch(fact<-ceiling(dt.ser/dt),
             error=fact<-1)
    #  if(fact != dt.ser/dt)
    # {
    #   warning("Time interval ", dt, "h is not a divisor of original interval ", dt.ser, "h")
    #  }
  }
  if(length(fact)>0)
  {
    if(fact>1)
    {
      # disaggregating observations from larger to smaller time interval by given factor
      # recalc interval
      dt <- dt.ser/fact
      #  times for new series
      #  tms <- seq(ser.start, ser.end, by=dt*3600)

      # new values duplicate series and rescale if not a rate in terms of a fixed
      # period e.g m/hr
      vals <- apply(as.matrix(ser), MARGIN=2, FUN=rep, each=fact)   # won't wotk if fact is not an integer
      vals <-matrix(coredata(vals), ncol=ncol(ser))

      tms <- seq(ser.start, along.with=vals, by=dt*3600)
      # if the value is a rate then it should be applied to all of the values in
      # the interval "as is". Otherwise each values needs to be divided across the
      # smaller time steps so that the total across the original time intervals is the same
      #if(!is.rate){ser.agg <- ser.agg/fact}
      ser.agg <- xts(vals, order.by=tms[1:nrow(vals)])

      names(ser.agg)<-names(ser)

      return(ser.agg)
    }
  }

  #  aggegrating to larger time interval and fewer observations
  return(aggregate_xts(ser, dt))
}

# aggregate given series to a lower time resolution
aggregate_xts <- function(ser,
                          dt,
                          fun=mean)  # use SUM if values are
{
  if(is.null(ser)){return(ser)}
  if(!is.zoo(ser)){stop("Time series input required")}
  ser.start <- start(ser)
  ser.end <- end(ser)
  #  qobs <- subset_zoo(ser, start, end)
  # frequency of the input series
  dt.ser <- diff(as.numeric(index(ser)))
  if(!all(dt.ser[] ==dt.ser[1]))
  {
    warning("Irregularly spaced time series supplied to aggregate_xts")
  }
  dt.ser <- modal(dt.ser)/3600  # in hours
  n.aggr <- dt/dt.ser
  if(length(n.aggr)>0)
  {
    if(n.aggr>1)
    {
      # aggregation (e.g quarter hour to hourly)
      tms <- seq(ser.start, ser.end, by=dt*3600)
      index.agg <- rep(tms, each=n.aggr)  # index(proj$rain)
      # trim - why needed ?
      index.agg <- index.agg[1:nrow(ser)]
      # aggregate using mean by default
      ser.agg <- aggregate(zoo(ser), by = index.agg, FUN=fun)
      names(ser.agg)<-names(ser)
      return(ser.agg)
    }
    else if(n.aggr<1)
    {
      fact <- 1/n.aggr
      return(disaggregate_xts(ser, dt=dt))
    }
  }
  return(ser)
}
