#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 FOR MONTHLY OPERATIONS:
#  countMonthlyRecords       Returns a series with monthly counts of records
#  rollMonthlyWindows        Returns start/end dates for rolling time windows
#  rollMonthlySeries         Rolls Monthly a 'timeSeries' on a given period
################################################################################


# DW:
# I think we should call these functions:
# countRecordsMonthly, rollWindowsMonthly, rollSeriesMonthly, ...


# ------------------------------------------------------------------------------


countMonthlyRecords <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a series with monthly counts of records

    # Example:
    #   x = as.timeSeries(data(msft.dat)); countMonthlyRecords(x)
    #   x = as.timeSeries(data(edhec)); countMonthlyRecords(x)

    # FUNCTION:
    
    # Check Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Check for Signal Series:
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # Count:
    ans <- rollMonthlySeries(x[, 1], period = "1m", by = "1m", FUN = NROW)
    colnames(ans) <- "Counts"
      
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation 

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rollMonthlyWindows <-
    function(x, period = "12m", by = "1m")
{   
    # A function implemented by  Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns start and end dates for rolling time windows

    # Arguments:
    #   x - a 'timeSerie's object of asset returns
    #   period - a character string denoting the length of the rolling
    #       window, e.g. "24m" means 24 months
    #   by - a character string denoting the shift of the rolling window,
    #       e.g. "3m" means one quarter

    # FUNCTION:
    
    # Check Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Check for Signal Series:
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # Get Window Parameter:
    periodLength <- as.numeric(substr(period, 1, nchar(period)-1))
    periodUnit <- substr(period, nchar(period), nchar(period))
    byLength <- as.numeric(substr(by, 1, nchar(by)-1))
    byUnit <- substr(by, nchar(by), nchar(by))
    stopifnot(periodUnit == "m")
    stopifnot(byUnit == "m")

    # Get Window Parameter:
    periodLength <- as.numeric(substr(period, 1, nchar(period)-1))
    periodUnit <- substr(period, nchar(period), nchar(period))
    byLength <- as.numeric(substr(by, 1, nchar(by)-1))
    byUnit <- substr(by, nchar(by), nchar(by))
    stopifnot(periodUnit == "m")
    stopifnot(byUnit == "m")

    # Make Windows - We expect monthly data records ...
    positions <- time(x)
    Positions <- unique(timeFirstDayInMonth(positions))
    numberOfPositions <- length(Positions)
    startDates <- Positions[1:(numberOfPositions-periodLength)]
    endDates <- Positions[(periodLength+1):numberOfPositions]-24*3600

    # Windows:
    windows <- list(from = startDates, to = endDates)
    attr(windows, "control") = c(start = start(positions), end = end(positions))

    # Return Value:
    windows
}


# ------------------------------------------------------------------------------


rollMonthlySeries <-
    function(x, period = "12m", by = "1m", FUN, ...)
{   
    # A function implemented by  Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Rolls monthly a 'timeSeries' on a given period

    # Arguments:
    #   x - a 'timeSerie's object of asset returns
    #   period - a character string denoting the length of the rolling
    #       window, e.g. "24m" means 24 months
    #   by - a character string denoting the shift of the rolling window,
    #       e.g. "3m" means one quarter
    #   FUN - function to be applied

    # FUNCTION:
    
    # Check Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Check for Signal Series:
    stopifnot(is(x, "timeSeries"))
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # Settings:
    windows <- rollMonthlyWindows(x = x[, 1], period = period, by = by)

    # Apply Function:
    ans <- applySeries(x = x, from = windows$from, to = windows$to,
        FUN = FUN, ...)
   
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
      
    # Return Value:
    ans
}


################################################################################

