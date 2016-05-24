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
#  show,timeSeries           Prints a 'timeSeries' object
#  print,timeSeries          Prints a 'timeSeries' object
#  .print.timeSeries         Called by function print,timeSerie
################################################################################


setMethod("show", "timeSeries",
    function(object)
    {
        # A function implemented by Diethelm Wuertz and Yohan Chalabi

        # Description:
        #   Print method for an S4 object of class "timeSeries"

        # FUNCTION:

        # Check records to get printed:
        maxRmetrics <- as.numeric(getRmetricsOptions("max.print"))
        maxR <- as.numeric(getOption("max.print"))
        maxR <- floor(maxR / (NCOL(object) + NCOL(object@recordIDs)))
        max <- min(na.omit(c(maxRmetrics, maxR, Inf)))
        #-> Inf to cast case when maxRmetrics and maxR are NULL

        if (ptest <- ((omitted <- NROW(object) - max) > 0))
            object <- object[seq.int(max),]

        data <- as(object, "matrix")
        recordIDs <- object@recordIDs
        FinCenter <- finCenter(object)

        # Series:
        cat(FinCenter, "\n", sep = "")
        if (prod(dim(recordIDs)) & (nrow(data) == NROW(recordIDs))) {
            dataIDs <- as.matrix(recordIDs)
            colnames(dataIDs) <- paste(colnames(dataIDs), "*", sep = "")
            #-> use format(data) to have same number of digits when timeSeries
            # is printed without @recordIDs
            print(cbind(format(data), dataIDs), quote = FALSE, right = TRUE)
        } else {
            print(data, quote = FALSE) #-> to be consistent with @recordIDs print
        }

        # print message
        if (ptest)
            cat(gettextf("...\n [ reached getRmetricsOption('max.print') | getOption('max.print') -- omitted %i rows ]]\n", omitted))

        # Return Value:
        invisible(NULL) # as specified in ?show
    }
)


# ------------------------------------------------------------------------------


.print.timeSeries <-
    function(x, FinCenter = NULL, format = NULL,
    style = c("tS", "h", "ts"), by = c("month", "quarter"), ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Allows for horizontal and ts like print output.

    # Arguments:
    #   x - an object of class timeSeries
    #   FinCenter - print with time stamps according to FinCenter
    #   format - use specified format for printing
    #   style - a character value specifying how to print:
    #       "tS" Rmetrics' default vertical print style
    #       "h" horizontal print style,
    #       "ts" R's base style for regular time series
    #   by - specifies the period for a regular time serie,
    #       note only active for style="ts".

    # Example:
    #   x = timeSeries(); print(x, format = "%d %b %Y"); x

    # FUNCTION:

    # Change Financial Center:
    if (!is.null(FinCenter))
        finCenter(x) <- FinCenter

    # Match Arguments:
    style = match.arg(style)
    by = match.arg(by)

    # Change Format:
    if (is.null(format)) {
        charvec = rownames(x)
    } else {
        ans = timeDate(charvec = rownames(x),
            zone = "GMT", FinCenter = finCenter(x))
        if (format == "%Q") {
            Quarters = rep(paste("Q", 1:4, sep = ""), each = 3)
            Y = atoms(ans)[, 1]
            Q = Quarters[atoms(ans)[, 2]]
            charvec = paste(Y, Q)
        } else {
            charvec = format(ans, format)
        }
    }

    # Styles:
    if (style == "tS") {
        cat(finCenter(x), "\n")
        X <- getDataPart(x)
        rownames(X) = charvec
        print(X)
    } else if (style == "h") {
        stopifnot(isUnivariate(x))
        # print(as.vector(x))
        ans = as.matrix(x)[,1]
        names(ans) = charvec
        print(ans)
    } else if (style == "ts") {
        freq = c(month = 12, quarter = 4)
        start(x)
        start = unlist(atoms(start(x)))
        end = unlist(atoms(end(x)))
        ts = ts(as.vector(x), start[1:2], end[1:2], freq[by])
        print(ts)
    }

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


setMethod("print", "timeSeries",
    .print.timeSeries)


################################################################################

