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
# FUNCTION:                 DESCRIPTION:
#  aggregate,timeSeries      Aggregates a 'timeSeries' object
# FUNCTION:                 DESCRIPTION:
#  daily2monthly             Aggregates a daily to monthly 'timeSeries' object
#  daily2weekly              Aggregates a daily to weekly 'timeSeries' object
################################################################################


.aggregate.timeSeries <- 
    function(x, by, FUN, ...)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #   Aggregates a 'timeSeries' object

    # Details:
    #   This function can be used to aggregate and coursen a
    #   'timeSeries' object.

    # Arguments:
    #   x - a 'timeSeries' object to be aggregated
    #   by - a calendarical block
    #   FUN - function to be applied, by default 'colMeans'
    #   ... - additional argument to be passed to the newly generated
    #       'timeSeries' object

    # Value:
    #   Returns a S4 object of class 'timeSeries'.

    # Examples:
    # Quarterly Aggregation:
    #   m = matrix(rep(1:12,2), ncol = 2)
    #   ts = timeSeries(m, timeCalendar())
    #   Y = getRmetricsOptions("currentYear"); Y
    #   from = paste(Y, "04-01", sep = "-"); to = paste(Y+1, "01-01", sep = "-")
    #   by = timeSequence(from, to, by = "quarter") - 24*3600; by
    #   ts; aggregate(ts, by, sum)
    # Weekly Aggregation:
    #   dates = timeSequence(from = "2009-01-01", to = "2009-02-01", by = "day")
    #   data = 10 * round(matrix(rnorm(2*length(dates)), ncol = 2), 1); data
    #   ts = timeSeries(data = data, charvec = dates)
    #   by = timeSequence(from = "2009-01-08",  to = "2009-02-01", by = "week")
    #   by = by - 24*3600; aggregate(ts, by, sum)

    # FUNCTION:

    # Check Arguments:
    if (!((inherits(by, "timeDate") && x@format != "counts") ||
          (is.numeric(by) && x@format == "counts")))
        stop("'by' should be of the same class as 'time(x)'", call.=FALSE)
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation

    # Make sure that x is sorted:
    if (is.unsorted(x))
        x <- sort(x)

    # Sort and remove double entries in by:
    by <- unique(sort(by))

    INDEX <- findInterval(x@positions, as.numeric(by, "sec") + 1)
    INDEX <- INDEX + 1
    is.na(INDEX) <- !(INDEX <= length(by))

    # YC : ncol important to avoid problems of dimension dropped by apply
    data <- matrix(apply(getDataPart(x), 2, tapply, INDEX, FUN), ncol=ncol(x))
    rownames(data) <- as.character(by[unique(na.omit(INDEX))])
    colnames(data) <- colnames(x)
    ans <- timeSeries(data, ...)

    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
      
    # Return Value:
    ans
    
}


setMethod("aggregate", "timeSeries", function(x, by, FUN, ...)
          .aggregate.timeSeries(x, by, FUN, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
aggregate.timeSeries <- function(x, ...) .aggregate.timeSeries(x, ...)


################################################################################


daily2monthly <-
  function (x, init = FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Converts daily to monthly series   
    
    # Arguments:
    #    x - daily time series
    #    init - should the index series converted to a wealth series
    
    # FUNCTION:
    
    # Save Colnames:
    colNames <- colnames(x)
    
    # Fill to end of Month:
    Time <- unique(sort(timeLastDayInMonth(time(x))))
    x.endOfMonth <- x[nrow(x), ]
    time(x.endOfMonth) <- rev(Time)[1]
    x <- rbind(x, x.endOfMonth)
    x <- alignDailySeries(x, include.weekends=TRUE)
    
    # Cut Properly on end of Month:
    today <- timeDate(Sys.Date())
    first <- timeFirstDayInMonth(today)
    x <- x[time(x) < first, ]
    Time <- unique(sort(timeLastDayInMonth(time(x))))
    
    # Align Properly:
    mSeries <- alignDailySeries(x, include.weekends=TRUE)
    mSeries <- mSeries[Time, ]
    
    # Optionally Initialize:
    if (init) 
        for (i in 1:ncol(mSeries)) mSeries[, i] <- mSeries[, 
            i]/as.vector(mSeries[1, i])
    colnames(mSeries) <- colNames
    
    # Return Value:
    mSeries
}


# -----------------------------------------------------------------------------


daily2weekly <- 
  function(x, startOn="Tue", init=FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Converts daily to weekly series   
    
    # Arguments:
    #    x - daily time series
    #    init - should the index series converted to a wealth series
    
    # FUNCTION:
    
    # Convert Series:
    mSeries <- alignDailySeries(x, include.weekends = TRUE)
    start <- which(dayOfWeek(time(mSeries[1:7, ])) == startOn)
    mSeries <- mSeries[seq(start, nrow(mSeries), by = 7), ]
    
    # Wealth Initialization:
    if (init) for (i in 1:ncol(mSeries)) 
        mSeries[, i] <- mSeries[, i]/as.vector(mSeries[1, i])
    
    # Return Value:
    mSeries
}


###############################################################################


