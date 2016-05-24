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


###############################################################################
# FUNCTION:                 DESCRIPTION:
#  diff,timeSeries           Differences a 'timeSeries' object
###############################################################################


.diff.timeSeries <-
  function(x, lag = 1, diff = 1, trim = FALSE, pad = NA, ...)
{
    # A function implemented by Diethelm Wuertz
    # Modified by Yohan Chalabi

    # Description:
    #   Differences 'timeSeries' objects.

    # Arguments:
    #   x - a 'timeSeries' object.
    #   lag - an integer indicating which lag to use.
    #       By default 1.
    #   diff - an integer indicating the order of the difference.
    #       By default 1.
    #   trim - a logical. Should NAs at the beginning of the
    #       series be removed?
    #   pad - a umeric value with which NAs should be replaced
    #       at the beginning of the series.

    # Value:
    #   Returns a differenced object of class 'timeSeries'.

    # FUNCTION:
  
    # Ceck Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation

    # Convert:
    y <- getDataPart(x) # as.matrix(x)

    # Check NAs:
    # if (any(is.na(y))) stop("NAs are not allowed in time series")

    # Difference:
    z <- diff(y, lag = lag, difference = diff)
    diffNums = dim(y)[1] - dim(z)[1]

    # Trim Positions:
    if (!trim) {
        zpad <- matrix(0*y[1:diffNums, ] + pad, nrow = diffNums)
        z <- rbind(zpad, z) }
    pos <-
        if (!trim)
            x@positions
        else
            x@positions[-(1:diffNums)]

    # Record IDs:
    df <- x@recordIDs
    if (trim && sum(dim(df)) > 0) {
        df <- df[-seq.int(diffNums), , drop = FALSE]
        rownames(df) <- seq.int(NROW(df))
    }

    # Diff Result:
    ans <- timeSeries(data = z, charvec = pos, units = colnames(z),
               format = x@format, zone = x@FinCenter,
               FinCenter = x@FinCenter, recordIDs = df)
  
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
  
    # Return Value:
    ans
}


# -----------------------------------------------------------------------------


setMethod("diff", "timeSeries",
          function(x, lag = 1, diff = 1, trim = FALSE, pad = NA, ...)
          .diff.timeSeries(x, lag, diff, trim, pad, ...)
          ##x <- getDataPart(x)
          ##callGeneric()
          )


# until UseMethod dispatches S4 methods in 'base' functions
diff.timeSeries <- function(x, ...) .diff.timeSeries(x, ...)


###############################################################################

