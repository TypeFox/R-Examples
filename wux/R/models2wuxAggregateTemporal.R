
# -----------------------------------------------------------------
# $Author: $
# $Date: $
# $Rev: $
# -----------------------------------------------------------------



## ----------------------------------------------------------------
## temporal aggregation functions:
##
## - AggregateTemporal
## - GetDayMonthYearFactor
## - GetMonthYearFactor
## - seq.date.daily
## - ConvertPeriodToFactor
## - tapply.over.array
## ----------------------------------------------------------------



AggregateTemporal <- function(data.array,
                              temporal.aggr,
                              nc.time.list,
                              start.date,
                              end.date,
                              ...) {

  ## aggregate 3-dim. data array (lon, lat, time) to temporal statistics
  ## specified in 'temporal.aggregation' in 'user.input'.
  ##
  ## Args:
  ##   data.array: 3-dimensional data array (lon, lat, time) to be aggregated.
  ##   temporal.aggr: list of temporal statistics specified in
  ##     temporal.aggregation in user.input.
  ##   nc.time.list: list of time information of the NedCDF files.
  ##   start.date: start date specified in user.input.  
  ##   end.date: end.date specified in user.input.
  ##   '...' : optional input. not in use right now.   
  ##
  ## Returns:
  ##   data.list: 3-dimensional data array (lon, lat, time) according to
  ##     temporal statistics specified in 'temporal.aggregation' in 'user.input'.
  ## 
  ## History:
  ##   2011-07-26 | original code (geh)
  ##   2014-04-10 | filling of missing values in the time series added (thm)

  ## get period factors of the data of the first statistics level
  time.factor.df <- GetMonthYearFactor(nc.time.list, start.date, end.date, ...)

  ## keep array dimensions
  array.dim <- dim(data.array)

  ## Here we need to fill the original data.array with NAs in case the
  ## timeseries is not complete 

  ## are we facing monthly or daily data?
  temporal.resolution <- nc.time.list[[1]][["time.step"]]
  if (! temporal.resolution %in% c("daily", "monthly"))
    stop("UNKNOWN TEMPORAL RESOLUTION OF FILE: ", temporal.resolution)
  ## and what kind of calendar?
  cal.in.netcdf <- nc.time.list[[1]]["calendar"]
  if (cal.in.netcdf %in% c("days", "gregorian", "julian", "empty", "standard", "no entry", "proleptic_gregorian"))
    calendar <- "julian"
  else if (cal.in.netcdf %in% c("360_day", "360_days", "360 days", "360days"))
    calendar <- "360"
  else if (cal.in.netcdf %in% c("365_day", "365_days", "365 days", "365days", "noleap", "no_leap"))
    calendar <- "365"
  else
    stop("UNKNOWN CALENDAR TYPE IN NETCDF FILE: ", cal.in.netcdf)
  
  ## get both date vectors from NetCDF file and the theoretical vector defined by user
  if (temporal.resolution == "daily"){
    ## get dates of NetCDF file (same dim as dim of array) - either daily or monthly
    date.as.is.fac <- GetDayMonthYearFactor(nc.time.list, start.date, end.date, ...)
    date.as.is <- paste(date.as.is.fac[,1], date.as.is.fac[,2], date.as.is.fac[,3], sep = "-" )
    ## get dates that SHOULD be used according to start.date and end.date declared in userinput
    ## by the user himself.
    dates <- seq.date.daily(start.date, end.date, calendar)
    year.as.should <- as.integer(format(dates, format = "%Y"))
    month.as.should <- as.integer(format(dates, format = "%m"))
    day.as.should <- as.integer(format(dates, format = "%d"))
    dates.as.should.be <- paste(year.as.should, month.as.should, day.as.should, sep = "-")
  } else if (temporal.resolution == "monthly") {
    ## get dates of NetCDF file (same dim as dim of array) - either daily or monthly
    date.as.is <- GetMonthYearFactor(nc.time.list, start.date, end.date, ...)
    date.as.is <- paste(date.as.is[,1], date.as.is[,2], sep = "-" )
    ## get dates that SHOULD be used according to start.date and end.date declared in userinput
    ## by the user himself. Calendar is not relevant anymore.
    dates <- seq(start.date, end.date, by = "month")
    year.as.should <- as.integer(format(dates, format = "%Y"))
    month.as.should <- as.integer(format(dates, format = "%m"))
    dates.as.should.be <- paste(year.as.should, month.as.should, sep = "-")
  }

  ## check if the date in NetCDF files was read out correctly, i.e. if
  ## there are any timeslices read in more than once
  if (length(unique(date.as.is)) != length(date.as.is))
    stop("THERE ARE DUPLICATE DATES IN YOUR NETCDF FILE. CHECK YOUR InitModelsDictionary AND DELETE DUPLICATE FILES. (or maybe you defined a wrong temporal resolution in your init.models dictionary config).")
  
  ## check if there are missing time slices (i.e. missing files) in
  ## the NetCDF file. 
  are.timesteps.missing <- !dates.as.should.be %in% date.as.is
  ## (HACK) in case of a 360 days calendar we do not perform this check
  ## on daily bases, merely on a monthly basis
  if (calendar == 360){
      dates.as.should.be.mon <- paste(year.as.should, month.as.should, sep = "-")
      date.as.is.mon <- paste(date.as.is.fac[,1], date.as.is.fac[,2], sep = "-" )
      are.timesteps.missing <- !dates.as.should.be.mon %in% date.as.is.mon
 }
      
 
  ## change dimension of data.array in case of missing values and fill with NAs
  if ( any(are.timesteps.missing) ){
    cat("WARNING: THERE ARE ", sum(are.timesteps.missing) ,
        " MISSING TIMESTEPS IN THE MODEL! THIS WILL PROBABLY YIELD NAs IN THE FINAL OUTPUT\n")
    ## fill missing timeslices with NA
    data.array.tmp <- array(NA, c(array.dim[1:2], length(are.timesteps.missing)))
    data.array.tmp[,,!are.timesteps.missing] <- data.array
    ## swap
    data.array <- NULL
    data.array <- data.array.tmp
    data.array.tmp <- NULL
    ## get new period factors of the data of the first statistics level 
    year.factor <- factor(format(dates, format = "%Y"))
    month.factor <- factor(format(dates, format = "%m"))
    time.factor.df <-
      data.frame(year.factor = year.factor,
                 period.factor = month.factor)
  } else {
    ## as usual... keep the array as it is
  }

  
  ## OK, now start the real agregation process...

  ## for loop over statistics levels
  for (stat.level in as.character(names(temporal.aggr))) {
    cat(paste("      level of temporal aggregation:",
              stat.level, "\n", ""))
    by.time <- temporal.aggr[[stat.level]][["period"]]

    ## convert desired period to factors
    time.factor.df[["period.factor"]] <-
      ConvertPeriodToFactor(as.character(time.factor.df[["period.factor"]]),
                            by.time)
    
    ## case of time series or climatological statistics
    if (temporal.aggr[[stat.level]][["time.series"]]) {
      ## error in case of year.factor is NULL
      if (is.null(time.factor.df$year.factor)) {
        cat(paste("ERROR: AGGREGATION OF TIME SERIES OF",
                  stat.level, "NOT POSSIBLE!"), "\n")
        cat(paste("CHECK user.input FOR", stat.level), "\n")
        stop()
      }

      ## in case of time series use combined factors of period and year
      use.factor <- factor(paste(time.factor.df[["year.factor"]],
                                 time.factor.df[["period.factor"]], ""))
      data.list <- tapply.over.array(data.array, use.factor, temporal.aggr[[stat.level]][["statistic"]], ...)

      ## reform list of data.array to array to call tapply.over.array
      ## in the next iteration
      data.array <- array(unlist(data.list), c(array.dim[1:2],
                                               length(names(data.list))))

      ## extract time information from generated list and generate new
      ## time factors
      time.list <- strsplit(names(data.list), " ")
      time.factor.df <- data.frame(year.factor = sapply(time.list, "[", 1),
                                   period.factor = sapply(time.list, "[", 2))

    } else {   
      ## in case of climatological statistics only use period factors
      use.factor <- time.factor.df[["period.factor"]]
      data.list <- tapply.over.array(data.array, use.factor,
                                     temporal.aggr[[stat.level]][["statistic"]], ...)

      ## reform list of data.array to array to call tapply.over.array
      ## in the next iteration
      data.array <- array(unlist(data.list), c(array.dim[1:2],
                                               length(names(data.list))))

      ## extract time information from generated list and generate new
      ## time factors
      time.factor.df <- data.frame(period.factor = names(data.list))
    }    
  }
  
  rm(data.array)
  return(data.list)
}



GetMonthYearFactor <- function(nc.time.list,
                               start.date,
                               end.date,
                               ...) {
  
  ## gets time vector of required NetCDF files and flags the
  ## time steps with corresponding year and period specified
  ## in 'temporal.aggregation' in 'user.input'
  ##
  ## Args:
  ##   nc.time.list: list of time information of the NedCDF files.
  ##   start.date: start date specified in user.input.
  ##   end.date: end.date specified in user.input.
  ##   filenames: Files from which to retrieve time-vector
  ##   '...': optional input. not in use right now.
  ##
  ## Returns:
  ##   time.factor.df: data frame with monthly and yearly factors
  ##     according to the NetCDF time vector
  ## 
  ## History:
  ##   2010-10-27 | original code (thm)
  ##   2011-07-26 | added function GetYearlyTimeSteps
  ##              | and storage of factors in a data frame (geh)

  ## extracting and concatenating required time and date vectors of NetCDF files
  time.vector.out <- NULL
  date.vector.out <- as.POSIXct(NA)
  nc.start.date <- nc.time.list[[1]][["start.date"]]
  calendar.type <-  nc.time.list[[1]][["calendar"]]
  for (ii in seq(along = nc.time.list)) {
    count <- nc.time.list[[ii]][["count"]]
    offset <- nc.time.list[[ii]][["offset"]]
    time.vector <- nc.time.list[[ii]][["time"]][offset:(offset + count - 1)]
    time.vector.out <- c(time.vector.out, time.vector)
    date.vector <- nc.time.list[[ii]][["date.time"]][offset:(offset + count - 1)]
    date.vector.out <- c(date.vector.out, date.vector)
  }
  date.vector.out <- date.vector.out[!is.na(date.vector.out)]

  ## flag time steps with corresponding months and years
  month.factor <- factor(as.integer(format(date.vector.out, format = "%m")))
  year.factor  <- factor(as.integer(format(date.vector.out, format = "%Y")))
 
  ## create and return data frame with monthly and yearly factors
  time.factor.df <- data.frame(year.factor = year.factor,
                               period.factor = month.factor)

  return(time.factor.df)
}

GetDayMonthYearFactor <- function(nc.time.list,
                                  start.date,
                                  end.date,
                                  ...) {
  
  ## gets time vector of required NetCDF files and flags the
  ## time steps with corresponding year and period specified
  ## in 'temporal.aggregation' in 'user.input'
  ##
  ## Args:
  ##   nc.time.list: list of time information of the NedCDF files.
  ##   start.date: start date specified in user.input.
  ##   end.date: end.date specified in user.input.
  ##   filenames: Files from which to retrieve time-vector
  ##   '...': optional input. not in use right now.
  ##
  ## Returns:
  ##   time.factor.df: data frame with monthly and yearly factors
  ##     according to the NetCDF time vector
  ## 
  ## History:
  ##   2014-04-09 | original code, copied from GetMonthYearFactor (thm) 

  ## extracting and concatenating required time and date vectors of NetCDF files
  time.vector.out <- NULL
  date.vector.out <- as.POSIXct(NA)
  nc.start.date <- nc.time.list[[1]][["start.date"]]
  calendar.type <-  nc.time.list[[1]][["calendar"]]
  for (ii in seq(along = nc.time.list)) {
    count <- nc.time.list[[ii]][["count"]]
    offset <- nc.time.list[[ii]][["offset"]]
    time.vector <- nc.time.list[[ii]][["time"]][offset:(offset + count - 1)]
    time.vector.out <- c(time.vector.out, time.vector)
    date.vector <- nc.time.list[[ii]][["date.time"]][offset:(offset + count - 1)]
    date.vector.out <- c(date.vector.out, date.vector)
  }
  date.vector.out <- date.vector.out[!is.na(date.vector.out)]

  ## flag time steps with corresponding months and years
  day.factor   <- factor(as.integer(format(date.vector.out, format = "%d")))
  month.factor <- factor(as.integer(format(date.vector.out, format = "%m")))
  year.factor  <- factor(as.integer(format(date.vector.out, format = "%Y")))
 
  ## create and return data frame with monthly and yearly factors
  time.factor.df <- data.frame(year.factor = year.factor,
                               period.factor = month.factor,
                               day.factor = day.factor)

  return(time.factor.df)
}



ConvertPeriodToFactor <- function(period.vector,
                                  by.time) {
  
  ## looks which  period.vector element containy 'by.time'
  ## element and assigns each element to the corresponding
  ## 'by.time' element.
  ##
  ## Args:
  ##   period.vector: numeric or character vector indicating the periods.
  ##   by.time: named list of numeric or character vectors containing the
  ##     aggregation period.
  ##
  ## Returns:
  ##   new.label.vector: vector of same length as 'period.vector' containing
  ##     the corresponding label from 'by.time'.
  ##
  ## History:
  ##   2010-10-27 | original code (thm)
  ##   2011-07-26 | some adaptations (geh)
  
  ## coerce to integer
  if (is.numeric(sapply(by.time, "[", 1))) {
    period.vector <- as.numeric(period.vector)
  }
  
  new.label.vector <- period.vector
  ## get one by one e.g. season 
  for (ii in seq(along = by.time)) {
    name <- names(by.time[ii])
    aggregate.period.vector <- by.time[[ii]]
    new.label.vector[new.label.vector %in% aggregate.period.vector] <- name
  }
  ## only make factor if there is a user.input 'by.time', else return 'NULL'
  if (length(by.time) == 0)
    new.label.vector <- NULL
  
  ## exit function
  return(new.label.vector)
}


  seq.date.daily <- function(startdate, enddate, calender.days="julian"){
    ## calculates daily date sequence from startdate to enddate, defined on
    ## specific calendar (calender.days).
    ##
    ## Args:
    ##   startdate:     POSIX time object.
    ##   enddate:       POSIX time object.
    ##   calender.days: Sring. which caledar should be used? options are
    ##                  "360" days, "365" days or "julian" (default).
    ##
    ## Returns:
    ##   Sequemce of daily POSIX date objects from startdate to enddate.
    ##   The length of the vector is defined by calendar.days.
    ## 
    ## History:
    ##   2014-04-10 | Original code (thm)

    ## start by complete timeseries on classical (julian) calendar
    inflated.seq <- seq(startdate, enddate, by="day")

    calender.days <- as.character(calender.days)

    ## which calendar?
    if (calender.days == "julian") {
      out.seq <- inflated.seq
    } else if (calender.days  %in% c("360", "365")){
      inflated.conv.seq <- conv365tojulian(seq(0, length(inflated.seq) - 1),
                                           startdate, calender.days=calender.days)
      out.seq <- inflated.conv.seq[enddate > inflated.conv.seq]
    } else {
      stop("UNKNOWN CALENDAR TYPE")
    }

    return(out.seq)
  }



tapply.over.array <- function(x, INDEX, fun, na.rm = FALSE, ...) {
  ## Helperfunction which applies 'tapply' on the c(1,2)-margin of a
  ## 3-dimensional array, taking the 'INDEX' vector as the ragged array
  ## for the time-dimension (i.e. the last dimension). See '?tapply' for more help. 
  ##
  ## Args:
  ##   x: 3-dimensional data array (lon-lat-time) to be aggregated.
  ##   INDEX: Factor with length of time-dimension.
  ##   fun: The function to be applied.
  ## 
  ## Returns:
  ##   List with 2D arrays with reduced time dimension according to INDEX-levels.
  ## 
  ## History:
  ##   2010-10-27 | Original code (thm)
  ##   2011-02-04 | na.rm=T added for mean calc (geh, thm)
  ##   2011-05-26 | keep dimension in lapply with drop=FALSE (thm)
  ##   2014-04-10 | na.rm behavior added (thm)
  ## 
  ## TODO(thm): how much memory does this function need? maybe
  ##   without z <- x[,,INDEX == y], and just plug in?
  
  stopifnot(length(INDEX) == dim(x)[3])

  ## aggregate
  aggr <- lapply(unique(INDEX), function(y) {z <- x[,,INDEX == y, drop = FALSE];
                                             apply(z, c(1,2), fun, na.rm = na.rm)})
  names(aggr) <- unique(INDEX)
  ## convert list to array (not needed any more, list are better to deal with)
  ##  aggr <- array(unlist(aggr), dim=c(dim(x)[1], dim(x)[2], length(aggr)))
  
  return(aggr)
}
