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
#  align                     Aligns a 'timeSeries object' to time stamps
#  .align.timeSeries         Aligns a 'timeSeries object' to time stamps
################################################################################


# DW: See also ...
#   in package timeDate
#   setMethod("align", "ANY",
#   setMethod("align", "timeDate",


# ------------------------------------------------------------------------------


.align.timeSeries <- 
  function(x, by = "1d", offset = "0s", 
    method = c("before", "after", "interp", "fillNA", "fmm", "periodic", 
      "natural", "monoH.FC"), include.weekends = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Aligns a 'timeSeries' object to equidistant time stamps

    # Arguments:
    #   x - an object of class "timeSeries".
    #   by - 
    #   offset -
    #   method -
    #       "before" - use the data from the row whose position is
    #           just before the unmatched position;
    #       "after" - use the data from the row whose position is
    #           just after the unmatched position;
    #       "linear" - interpolate linearly between "before" and
    #           "after".
    #       "fillNA" - fill missing days with NA values
    #   include.weekends - a logical value. Should weekend dates be
    #       included or removed?

    # Example:
    #      data(usdthb)
    #      data = matrix(usdthb[, "BID"])
    #      charvec = as.character(usdthb[, "XDATE"])
    #      USDTHB = timeSeries(data, charvec, format = "%Y%M%d%H%M")
    #      align(USDTHB, by = "3h", offset = "92m")
    #      MSFT = as.timeSeries(data(msft.dat))
    #      align(MSFT)

    # See also:
    #   in package timeDate
    #   setMethod("align", "ANY",
    #   setMethod("align", "timeDate",

    # FUNCTION:
    
    # Settings:
    Title <- x@title
    Documentation <- x@documentation

    # Check for Signal Series:
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # check if series sorted
    if (is.unsorted(x)) x <- sort(x)

    # Adjustment:
    Method <- match.arg(method)
    fun <- switch(Method,

                  before = function(x, u, v)
                  approxfun(x = u, y = v, method = "constant", f = 0, ...)(x),

                  after = function(x, u, v)
                  approxfun(x = u, y = v, method = "constant", f = 1, ...)(x),

                  interp = ,
                  fillNA = function(x, u, v)
                  approxfun(x = u, y = v, method = "linear", f = 0.5, ...)(x),

                  fmm = ,
                  periodic = ,
                  natural = ,
                  monoH.FC = function(x, u, v)
                  splinefun(x = u, y = v, method = Method, ...)(x))

    # Align timeDate stamps
    td1 <- time(x)
    td2 <- align(td1, by = by, offset = offset)

    # extract numerical values
    u <- as.numeric(td1, units = "secs")
    xout <- as.numeric(td2, units = "secs")

    N = NCOL(x)
    data <- matrix(ncol = N, nrow = length(td2))
    xx <- getDataPart(x)
    for (i in 1:N) {

        v <- as.vector(xx[, i])

        # New Positions and approximated values:
        yout <- fun(xout, u, v)
        if (Method == "fillNA") yout[!(xout %in% u)] = NA

        # Compose data:
        data[, i] = yout
    }

    # build time series
    ans <- timeSeries(data, td2, units = colnames(x))

    # Remove Weekends:
    if(!include.weekends) ans <- ans[isWeekday(td2), ]
    
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation

    # Return Value:
    ans
}

# ------------------------------------------------------------------------------


setMethod("align", "timeSeries", .align.timeSeries)


################################################################################

