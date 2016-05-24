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
#  lag,timeSeries            Lags a 'timeSeries' object
################################################################################


setMethod("lag" , "timeSeries",
    function(x, k = 1, trim = FALSE, units = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Lags a 'timeSeries' object.

    # Arguments:
    #   x - a 'timeSeries' object.
    #   k - an integer indicating which lag to use. By default 1.
    #       Note, negative lags are to data in the future.
    #   trim - a logical. Should NAs at the beginning of the
    #       series be removed? By default FALSE.
    #   units -
    #   ... -

    # Details:
    #   The arguments differ in the following way from the function
    #   stats::lag - lag(x, k = 1, ...)

    # Value:
    #   Returns a lagged object of class 'timeSeries'.

    # Example:
    #   SPI = 100* as.timeSeries(data(LPP2005REC))[1:20, "SPI"]
    #   Note negative lags are to data in the future !
    #   lag(SPI, k = -2:2)
    #   lag(SPI, k = 0:2 , trim = TRUE)

    # FUNCTION:

    # Check Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Internal Function:
    tslagMat <- function(x, k = 1) {
        # Internal Function:
        tslag1 = function(x, k) {
            y = x
            if (k > 0) y = c(rep(NA, times = k), x[1:(length(x)-k)])
            if (k < 0) y = c(x[(-k+1):length(x)], rep(NA, times = -k))
            y }
        # Bind:
        ans <- NULL
        for (i in k) {
            ans <- cbind(ans, tslag1(x, i)) }
         # As Vector:
        if (length(k) == 1) ans <- as.vector(ans)
        # Return Value:
        ans }

    # Convert:
    y <- getDataPart(x)
    Dim <- dim(y)[2]

    # Lag on each Column:
    z <- NULL
    for (i in 1:Dim) {
         ts <- tslagMat( y[, i], k = k)     #, trim = FALSE)
         z <- cbind(z, ts) }

    # Positions
    pos <- x@positions

    # Record IDs:
    df <- x@recordIDs

    # Trim:
    if (trim){
        idx <- !is.na(apply(z, 1, sum))
        z <- z[idx, , drop = FALSE]
        pos <- pos[idx]
        if (sum(dim(df)) > 0) {
            df <- df[idx, , drop = FALSE]
            rownames(df) <- seq.int(sum(idx))
        }
    }

    # Augment Colnames:
    cn <- colnames(x)
    a <-
        if (is.null(units))
            # ensure that colnames is replicated according to the length
            # of lag indexes.
            as.vector(matrix(cn, nrow = length(k), ncol = length(cn), byrow = TRUE))
        else
            units

    kcols <- rep(k, times = ncol(y))
    b <- paste("[", kcols, "]", sep="")
    ab <- paste(a, b, sep = "")
    units <- ab
      
    # Result:
    ans <- timeSeries(data = z, charvec = pos, units = units,
        format = x@format, FinCenter = x@FinCenter, recordIDs = df)
     
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
      
    # Return Value:
    ans

})


# until UseMethod dispatches S4 methods in 'base' functions
lag.timeSeries <- function(x, ...) timeSeries::lag(x, ...)


################################################################################

