# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2014-04-11 10:24:47 +0200 (Fri, 11 Apr 2014) $
# $Rev: 273 $
# ----------------------------------------------------------------

## ----------------------------------------------------------------
## You can find tools for manipulating date- and time vectors here:
##
## - Read time variable from NetCDF and ...
## - ... extract its month components.
## - Labeling months with summer/winter, or seasons.
## ----------------------------------------------------------------

ReadNetCdfTimeData <- function(filenames,
                               what.timesteps = NULL,
                               calendar = NULL,
                               time.units = NULL,
                               count.first.time.value = NULL,
                               ...) {
  ## Reads time vector from 'filenames', converts it to POSIXt date vector
  ## and stores them in lists.
  ## Note: When dealing with daily data, 360 days and 365 days calendars will
  ## return a POSIXt date vector with missing 5,6 (with 360days) or 0/1  
  ## (365days) days at the end of every year. This is not really a restriction,
  ## but in future we should think of a 365 days and 360 days POSIX class.
  ##
  ## Args:
  ##   filenames: Character vector of filenames for which to read time array.
  ##   calendar: Character indicating calendar used in climate simulation.
  ##             Right now possible: "julian" for usual calender or "360_days".
  ##             If set NULL (default), this function tries to retrieve the
  ##             calender attribute from the NetCDF file.
  ##   time.units: Character indicating time units stored in NetCDF file.
  ##             Right now possible: "days since ...." or "seconds since ...".
  ##             If set NULL (default), this function tries to retrieve the
  ##             time.units attribute from the NetCDF file.
  ##   count.first.time.value: In case the file calender has 360 days but the
  ##                           first timestep from the NetCDF time-variable has
  ##                           to be added in a "julian way", this argument must
  ##                           be flagged "julian". Else nothing happens.
  ##   '...': No additional meaning right now.
  ##
  ## Returns:
  ##   Named list of NetCDF files with time variable stored as (I) timesteps
  ##   (as in the NetCDF file) and (II) as POSIXct vector for smooth usage in R.
  ##
  ## Example:
  ##   infile <- "/data/reloclim/rcm/ENSEMBLES/MPI/REMO/ECHAM5-r3_A1B/
  ##              direct/MPI-M-REMO_SCN_ECHAM5_DM_25km_1961-1970_tas.nc"
  ##   nc.time <- ReadNetCdfTimeData(infile)
  ##   str(nc.time)
  ##
  ## History:
  ##   2010-10-27 | Original code (thm)
  ##   2010-11-25 | Keywords "calendar" and "time.units" added (thm)
  ##   2010-11-29 | Keyword "count.first.time.value" added (thm)
  ##   2012-03-14 | time in "minutes from...." added (seb)
  ##   2016-01-13 | change to "ncdf4" library

  cat("    reading time variable from NetCDF \n")
  ## generate empty time-series
  if ( is.list(filenames) ) {
    period.names <- as.vector(sapply(filenames, "[", 1))
    filenames <- names(filenames)
    write.period = TRUE
  } else {
    write.period = FALSE
  }
  list.return <- vector(mode = "list",  length(filenames))
  names(list.return) <- filenames

  ## reads one file after another and then concatenates the timeseries
  for ( file.index in c(1:length(filenames)) ) {

    list.return[[filenames[file.index]]] <-
      list(time = NULL, date.time = NULL, units = NULL,
           start.date = NULL, calendar = NULL,
           intervall.begin = NULL, intervall.end = NULL,
           count = NULL, offset = NULL)

### Reading data
    ## Start NetCDF reading
    ## nc <- open.ncdf(filenames[file.index]) old...
    nc <- nc_open(filenames[file.index])
     ## read dimensions of netcdf file
    ##nc.dims <- GetNetcdfDims(nc)
    ## read time unit-attribute in file
    ## units.attr <- att.get.ncdf(nc, "time", "units")
    units.attr <- ncatt_get(nc, "time", "units")
    ## read time steps
    ## TODO(thm): nicht hartcodiert sondern dynamisch nach vals suchen!!!
    ## timesteps <- get.var.ncdf(nc, "time")
    timesteps <- ncvar_get(nc, "time")
    ## flag if numbers of days in year is 360 or regular calendar
    ## calendar.attr      <- att.get.ncdf(nc, "time", "calendar")
    calendar.attr      <- ncatt_get(nc, "time", "calendar")
    have.calendar.attr <- calendar.attr$hasatt
    ## Close NetCDF file
    ## close.ncdf(nc)
    nc_close(nc)

### Getting data attributes
    ## getting time units
    if (!is.null(time.units)) {
      ## case we place time units by hand
      units <- time.units
    } else {
      have.units.attr <- units.attr$hasatt
      if (have.units.attr) {
        units <-  units.attr$value
      } else {
        stop("NETCDF: VARIABLE TIME HAS NO UNITS ATTRIBUTE. THATS BAD. ")
      }
    }
    ## getting calendar in case we didn't pass it through the arguments manually
    if (is.null(calendar)) {
      if (have.calendar.attr) {
        calendar <- calendar.attr$value
      } else {
        calendar <- "no entry"
      }
    }

    ## Time is: days since...
    ## here stringseperate... "days since"
    if (pmatch("days since", units, nomatch=0)) {
      start.date <- as.POSIXct(strsplit(units, "days since ")[[1]][2],
                               tz="GMT")
    }

    ## Time is: hours since
    if (pmatch("hours since", units, nomatch=0)) {
      start.date <- as.POSIXct(strsplit(units, "hours since ")[[1]][2],
                               tz="GMT")
      ## transform to daily basis (to be transformed back to hours...stupid?)
      ## TODO(thm): do not calculate on daily basis, but on ... second-basis?
      ## depend on what is best...
      timesteps <- timesteps / 24
    }
    
    ## Time is: seconds since...
    if (pmatch("seconds since", units, nomatch=0)) {
      start.date <- as.POSIXct(strsplit(units, "seconds since ")[[1]][2],
                               tz="GMT")
      ## transform to daily basis (to be transformed back to seconds...stupid?)
      ## TODO(thm): do not calculate on daily basis, but on ... second-basis?
      ## depend on what is best...
      timesteps <- timesteps / 60 / 60 / 24
    }
    ##Time is: minutes since
     if (pmatch("minutes since", units, nomatch=0)) {
      start.date <- as.POSIXct(strsplit(units, "minutes since ")[[1]][2],
                               tz="GMT")
     
      timesteps <- timesteps /  60 / 24
    }
    ## Time is: months since...
    if (pmatch("months since", units, nomatch=0)) {
      start.date <- as.POSIXct(strsplit(units, "months since ")[[1]][2],
                               tz="GMT")
    }

    
    ## critical case: file has 360 days calender BUT the first timestep adds
    ## to the time units like the file had julian date or noleap-date.
    ## eg:
    ##   file starts from 2031-01-01 and has 360 days calender
    ##   time:units = "days since 1950-01-01 00:00:00" ;
    ##   time = 29585.5, 29586.5, 29587.5, 29588.5, ...;
    ## 
    ## solution: define new startdate (julian adding) and adapt the
    ##           timesteps-vector

    if (pmatch("days since", units, nomatch=0) |
        pmatch("hours since", units, nomatch=0) |
        pmatch("seconds since", units, nomatch=0)|
        pmatch( "minutes since",units, nomatch=0 )) {
      if (calendar %in% c("360_day", "360_days", "360 days", "360days")) {
        if (!is.null(count.first.time.value))
          if (count.first.time.value == "julian") {
            start.date <- start.date + timesteps[1] * 60 * 60 * 24
            timesteps <- timesteps - timesteps[1]
          } else {
            if (count.first.time.value == "noleap") {
              start.date <- conv365tojulian(timesteps[1], start.date)
            }
          }
        date.time <- conv365tojulian(timesteps, start.date, calender.days="360")
      }

      ## 365 days calendar
      if (calendar %in% c("365_day", "365_days", "365 days", "365days", "noleap"
                          , "no_leap")) {
        if (!is.null(count.first.time.value))
          stop("count.first.time.value NOT DEFINED FOR noleap calendar.
FEEL FREE TO IMPLEMENT IT")

        date.time <- conv365tojulian(timesteps, start.date)
      } else {
        if (calendar %in% c("days", "gregorian", "julian", "empty", "standard", "no entry", "proleptic_gregorian")) {
          if (!is.null(count.first.time.value))
            stop("count.first.time.value NOT DEFINED FOR gregorian calendar.
FEEL FREE TO IMPLEMENT IT")

          ## normal gregorian calendar
          date.time <- start.date + timesteps * 60*60*24
        }
      }
    }

    if (pmatch("months since", units, nomatch=0)) {
      date.time.begin <-
        seq(start.date, by = "month", length = floor(timesteps[1] + 1))
      date.time <- seq(date.time.begin[length(date.time.begin)],
                       by = "month", length = length(timesteps))
    }

    ## Case: Time is specific format (as e.g. "day as %Y%m%d.%f") and
    ## not e.g. "days since 1949-12-01". Date time vector is simply
    ## the time vector reformated.
    ## cat("HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK")
    if (pmatch("day as", units, nomatch=0)) {
       units.format <- strsplit(units, "day as ")[[1]][2]
       if (length(grep("%.f", units.format)) == 0)
         ## if I find "%.f" in time units I'll delete it, as I don't
         ## know what it means and R can't interpret it.
         units.format <- sub('\\.\\%f', '', units.format)
       date.time <- as.POSIXct(strptime(timesteps, format = units.format, tz = "GMT"))
       start.date <- date.time[1]
     }

    if (is.null(date.time))
      stop("CALENDER NOT FOUND")

    ## check if time step is daily or monthly
    if( is.null(what.timesteps) ) {
      time.step <- "daily"
    } else {
      time.step <- what.timesteps
    }

    ## concatenate time vector
    list.return[[filenames[file.index]]][["time"]] <- timesteps
    list.return[[filenames[file.index]]][["date.time"]] <- date.time
    list.return[[filenames[file.index]]][["calendar"]] <- calendar
    list.return[[filenames[file.index]]][["units"]] <- units
    list.return[[filenames[file.index]]][["start.date"]] <- start.date
    list.return[[filenames[file.index]]][["time.step"]] <- time.step
    if ( write.period )
      list.return[[filenames[file.index]]][["period"]] <- period.names[file.index]
  }

  return(list.return)
}


GetTimeCountsOffsets <- function(nc.time, startdate, stopdate, check.period) {
  ## Calculates the NetCDF time offset and count from 'startdate' to 'stopdate'.
  ## It updates the list returned from ReadNetCdfTimeData by count and offset.
  ##
  ## Args:
  ##   nc.time: Object returned by ReadNetCdfTimeData. This is a list of NetCDF
  ##            filenames with the corresponding time variable stored.
  ##   startdate: POSIX time object. First date which should be considered for
  ##              offset. Note, that the time vector in 'nc.time' is stored as
  ##              POSIXct YYYY-MM-DD 12:00:00, thus days are considered being
  ##              at 12AM noon. If you want your offset start with a particular
  ##              day, make sure you are < 12:00:00 that day (see example below)
  ##   stopdate:  POSIX time object. Last date which should be considered for
  ##              count. Note, that the time vector in 'nc.time' is stored as
  ##              POSIXct YYYY-MM-DD 12:00:00, thus days are considered being
  ##              at 12AM noon. If you want your count end on a particular
  ##              day, make sure you are > 12:00:00 that day (see example below)
  ##  check.period: indicates the actual period (character e.g. 2021-2050)
  ##
  ## Returns:
  ##   Updated object (which has previously been retuned by ReadNetCdfTimeData)
  ##   with count and offset (for time dimension).
  ##
  ## Example:
  ## infile <- "/data/reloclim/rcm/ENSEMBLES/MPI/REMO/ECHAM5-r3_A1B/direct/MPI-M-REMO_SCN_ECHAM5_DM_25km_1961-1970_tas.nc"
  ## nc.time <- ReadNetCdfTimeData(infile)
  ## str(nc.time)
  ## startdate <- strptime("1961-01-01 05:00:00",format = "%Y-%m-%d %H:%M:%S")
  ## stopdate <- strptime("1961-01-23 23:00:00", format = "%Y-%m-%d %H:%M:%S")
  ## nc.time <- GetTimeCountsOffsets(nc.time, startdate, stopdate)
  ## str(nc.time)
  ##
  ## History:
  ##  2010-10-27 | Original code (thm)
  ##  2013-10-08 | extended overlapping time check to daily checking (thm)
  ##  2014-04-03 | extended overlapping time check to daily checking
  ##               erased as it made trouble for monthly data (thm)

  list.time.ret <- nc.time

  for (ii in (1:length(names(nc.time)))) {
    file.name <- names(nc.time)[ii]

    ## get count and offset
    date.seq <- list.time.ret[[file.name]][["date.time"]]
    first <- startdate <= date.seq
    last <- stopdate >= date.seq

    if (!(any(first))) {
      first.index <- NULL
    } else {
      first.index <- min(which(first == TRUE))
    }
    if (!(any(last))) {
      last.index <- NULL
    } else {
      last.index <- max(which(last == TRUE))
    }

    if (!is.null(first.index))
      offset <- first.index
    if (is.null(last.index)) {
      count <- NA
    } else {
      count <- last.index - first.index + 1
    }

    ## if file doesn't contain start- and enddate, set count and offset to NA
    first.datestep <- date.seq[1]
    last.datestep <- tail(date.seq, 1)
    if (startdate > last.datestep | stopdate < first.datestep)
      count <- offset <- NA

    ## store start- and stopdate from NetCDF file to "intervall.begin" and
    ## "intervall.end" tag
    if (is.na(offset)) {
      list.time.ret[[file.name]][["intervall.begin"]] <- NA
    } else {
      list.time.ret[[file.name]][["intervall.begin"]] <- date.seq[offset]
    }
    if (is.na(offset + count)) {
      list.time.ret[[file.name]][["intervall.end"]] <- NA
    } else {
      list.time.ret[[file.name]][["intervall.end"]] <-
        date.seq[offset + count - 1]
    }
    list.time.ret[[file.name]][["offset"]] <- offset
    list.time.ret[[file.name]][["count"]] <- count
    ## check for overlapping time periods of two consecutive files
    ## and correct the "count" and "interval.end" of the first file
    ## if necessary
    if (ii > 1) {
      if ( !is.null(check.period) ) {
        format.startdate <- as.character(format(startdate, format="%Y-%m"))
        check.startdate.first <-
          format.startdate %in% format(list.time.ret[[ii-1]]$date.time,
                                       format="%Y-%m")
        check.startdate.second <-
          format.startdate %in% format(list.time.ret[[ii]]$date.time,
                                       format="%Y-%m")
        if ( check.startdate.first & check.startdate.second ) {
          which.period <- c(list.time.ret[[ii-1]]$period,
                            list.time.ret[[ii]]$period) == check.period
          if ( which.period[1] ) {
            list.time.ret[[ii]]$count <- NA
            list.time.ret[[ii]]$offset <- NA
          }
          if ( which.period[2] ) {
            list.time.ret[[ii-1]]$count <- NA
            list.time.ret[[ii-1]]$offset <- NA
          }
        }
      }

      ## check if PREVIOUS file dates overlap with current file
      ## dates. This does not check if other files before the previous
      ## one overlap! In this case an error will be thrown in the
      ## function "AggregateTemporal" as it will detect duplicate
      ## dates. The user has to take care of this by himself in
      ## InitmodelsDict file and delete duplicates.
      if ( !is.na(list.time.ret[[ii-1]]$count) &
           !is.na(list.time.ret[[ii-1]]$offset) &
           !is.na(list.time.ret[[ii]]$count) &
           !is.na(list.time.ret[[ii]]$offset) ) {

        last.date <- format(list.time.ret[[ii]]$date.time[1], format="%Y-%m")
        is.overlap.vector <-
          last.date  == format(list.time.ret[[ii - 1]]$date.time, format="%Y-%m")

        if( !is.na(which(is.overlap.vector == TRUE)[1]) &
           last.date < format(stopdate, format="%Y-%m") ) {
          ## warning message
          cat("\n    !!!!! WARNING: OVERLAPPING TIME PERIODS IN NetCDF files. OVERWRITING",
                    "COUNT OF FILE: !!!!!", "\n", "")
          cat(paste("    ", names(nc.time)[ii - 1], "\n\n", sep=""))

          ## modify "count"
          list.time.ret[[ii - 1]]$count <- which(is.overlap.vector == TRUE) -
            list.time.ret[[ii - 1]]$offset
          ## modify "intervall.end"
          list.time.ret[[ii - 1]]$intervall.end <-
            list.time.ret[[ii - 1]]$date.time[which(is.overlap.vector == TRUE)
                                              - 1]

          ## check if "intervall.end" date is earlier than "intervall.begin"
          ## date
          diff.time <- as.numeric(list.time.ret[[ii - 1]]$intervall.end) -
            as.numeric(list.time.ret[[ii - 1]]$intervall.begin)
          if (diff.time < 0) {
            list.time.ret <- list.time.ret[ii]
          }
        }
      }

    }

  }

  return(list.time.ret)
}



conv365tojulian <- function(time, since, calender.days="365") {
  ## Converts timesteps from NetCDF file having 365 or 360 days calendar to
  ## POSIXt class timevector. The reason for this function is, that cannot simply
  ## add the timevector to the POSIX date we wish, as this would lead to complete
  ## nonsense data! Therefore we add the years and days seperately and leave out
  ##  may 31st, july 31st aug 31st, oct 31st and dec 31st (360 days) and
  ## march 31st if leapyear.
  ## For 365 days calendar we leave out feb 29th if being leap year.
  ##
  ## Args:
  ##   time: A vector with timesteps, which once was the time variable of a
  ##         NetCDF file. These timesteps have to be on DAILY basis!
  ##   since: POSIXt object, time starts with.
  ##   calender.days: Character or numeric. Calendar type. Can be "365" or "360"
  ##
  ## Returns:
  ##   POSIXct vector of dates within NetCDF file, having same length as NetCDF
  ##   time dimension.
  ##
  ## History:
  ##   2010-11-29 | Original code (thm)
  ##   2011-12-01 | Deleted hack, which was just not correct (it shifted one day back
  ##                ward, which for monthly data didnt matter, so data are processed
  ##                the same way) (thm)
  ##   2011-12-16 | handle properly the cases of 360 days and 365 days (thm)
  ##   2014-04-09 | changed to modulo calculation due to rounding error (thm)

  IsLeapYear <- function(year){
    ## small helperfunction returning boolean vector leapyear true/false
    ## http://en.wikipedia.org/wiki/Leap_year
    return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
  }

  ## convert arguments to useful format
  calender.days <- as.numeric(calender.days)
  since <- as.POSIXlt(since)

  ## get vector containing the "years since" and "day in that year"
  n.of.years <- time/calender.days
  year <- floor(n.of.years)

  ## merge "years since" and "day in that year" to a POSIXct date vector
  date.posix <- since
  ## get years since ... (by adding years)
  date.posix$year <- since$year + year
  ## boolean vector indicating if year is leapyear
  is.leapyear <- IsLeapYear(as.integer(format(date.posix, format = "%Y")))

  ## getting days within this year
  days <- time %% calender.days  ## OLD calculation: (n.of.years - year) * calender.days
  ## now comes the crucial part:
  ## we shift the days for 360days calendar and 365days calendar so, that:
  ## case 360days: all 31st (except jan and mar) disappear (for leapyears all
  ##               except jan)
  ## case 365days: for leapyears feb 29th disappears
  if (calender.days == 360){
    days[days >= 150 & !is.leapyear] <- days[days >= 150 & !is.leapyear] + 1#may 31th
    days[days >= 211 & !is.leapyear] <- days[days >= 211 & !is.leapyear] + 1#jul 31th
    days[days >= 242 & !is.leapyear] <- days[days >= 242 & !is.leapyear] + 1#aug 31th
    days[days >= 303 & !is.leapyear] <- days[days >= 303 & !is.leapyear] + 1#oct 31th
    days[days >= 364 & !is.leapyear] <- days[days >= 364 & !is.leapyear] + 1#dec 31th

    days[days >= 90  & is.leapyear] <- days[days >= 90  & is.leapyear] + 1 #mar 31th
    days[days >= 151 & is.leapyear] <- days[days >= 151 & is.leapyear] + 1 #may 31th
    days[days >= 212 & is.leapyear] <- days[days >= 212 & is.leapyear] + 1 #jul 31th
    days[days >= 243 & is.leapyear] <- days[days >= 243 & is.leapyear] + 1 #aug 31th
    days[days >= 304 & is.leapyear] <- days[days >= 304 & is.leapyear] + 1 #oct 31th
    days[days >= 365 & is.leapyear] <- days[days >= 365 & is.leapyear] + 1 #dec 31th
  }
  else if (calender.days == 365){
    days[days >= 60  & is.leapyear] <- days[days >= 60  & is.leapyear] + 1 #feb 29th
  }
  else {
    stop("INVALID CALENDAR TYPE IN conv365tojulian")
  }

  ## add days to POSIX object
  date.posix <- date.posix + days * 24 * 3600

  ## short output summary
  n.years.approx <- as.integer(strftime( tail(date.posix, 1), "%Y")) -
    as.integer(strftime( date.posix[1], "%Y")) + 1
  ## cat("Range in NetCDF file: ", as.character(date.posix[1]), " to ",
  ##     as.character(tail(date.posix, 1)),  " (~",
  ##     n.years.approx, " years)\n", sep="")

  return(date.posix)
}
