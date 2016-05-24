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
#  na.omit,timeSeries        Handles missing values in objects
#  .naOmitMatrix             Internal function called from na.omit,timeSeries
# OLD FUNCTIONS:            DESCRIPTION:
#  removeNA                  Remove NAs from a matrix object
#  substituteNA              Substitute NAs by zero, the column mean or median
#  interpNA                  Interpolate NAs using R's "approx" function
################################################################################


# DW:
# I think we should deprecate the following functions:
# removeNA, substituteNA, and interpNa since the function 
# na.omit() can already handle all these cases.


# DW:
# note we do interpolation with approx(), zoo also offers
# interpolation by splines, we should also add this.


################################################################################


.na.omit.timeSeries <- 
    function(object, method = c("r", "s", "z", "ir", "iz", "ie"),
        interp = c("before", "linear", "after"), ...)
{
    # Description
    #   Handles NAs in timeSeries objects

    # Detials:
    #   Linear Interpolation is done by the function approx.
    #   Spline interpolation like in zoo is not yet supported.

    # Arguments:
    #   object - an object of class 'timeSeries'
    #   method - how to handle NAs
    #   interp - how to interpolate NAs
    #   ... - arguments passed to function approx()

    # FUNCTION:

    # Check Arguments:
    stopifnot(is.timeSeries(object))
      
    # Extract Title and Documentation:
    Title <- object@title
    Documentation <- object@documentation
      
    # Settings:
    method <- match.arg(method)
    interp <- match.arg(interp)

    # Skip ?
    if (method == "s") return(object)

    # Handle NAs:
    if (method == "r") {
        # Remove NAs:
        # DW:
        # object <- stats:::na.omit.default(object)
        object <- as.timeSeries(na.omit(series(object)))
    } else if (method == "z") {
        # Substitute NAs by Zero's:
        object[is.na(object)] <- 0
    } else if (substr(method, 1, 1) == "i") {
        # Interpolate:
        interp = match.arg(interp)
        f = 0
        if (interp == "before") {
            interp = "constant"
            f = 0
        }
        if (interp == "after") {
            interp = "constant"
            f = 1
        }
        n = nrow(object)
        for (i in 1:ncol(object)) {
            y = object[, i]
            idy = (1:n)[!is.na(y)]
            ## DW: ... added
            # Linear/Constant Interpolation:
            y = approx(x = idy, y = y[idy], xout = 1:n, method = interp,
                f = f, ...)$y
            object[, i] = y
        }
        modID = FALSE
        if (method == "ir") {
            # Remove Start and End NAs:
            # DW:
            # object <- stats:::na.omit.default(object)
            object <- as.timeSeries(na.omit(series(object)))
        } else if (method == "iz") {
            # Set Start and End NAs to Zero:
            object[is.na(object)] = 0
        } else if (method == "ie") {
            n = nrow(object)
            for (k in 1:ncol(object)) {
                y  = object[, k]
                if (NA %in% y) {
                    start = sum(cumprod(is.na(y)))
                    if (start > 0) for (i in start:1) y[i] = y[i+1]
                    end = n+1 - sum(cumprod(rev(is.na(y))))
                    if (end <= n) for (i in end:n) y[i] = y[i-1]
                    object[, k] = y
                }
            }
        }
    }

    # Handle recordIDs ...
    recordIDs <- object@recordIDs
    modID <- c(r = TRUE, z = FALSE, ir = TRUE, iz = FALSE, ie = FALSE)
    if(modID[method] > 0 && sum(dim(recordIDs)) > 0 ) {
        index <- attr(object, "n.action")
        recordIDs <- recordIDs[index, ]
    }
      
    # Preserve Title and Documentation:
    object@title <- Title
    object@documentation <- Documentation 

    # Return Value:
    object
}


setMethod("na.omit", "timeSeries", function(object, 
    method = c("r", "s", "z", "ir", "iz", "ie"),
    interp = c("before", "linear", "after"), ...)
.na.omit.timeSeries(object, method, interp, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
na.omit.timeSeries <- function(object, ...) .na.omit.timeSeries(object, ...)


# ------------------------------------------------------------------------------


.naOmitMatrix <- 
    function(object, method = c("r", "s", "z", "ir", "iz", "ie"),
        interp = c("before", "linear", "after"))
{
    # Description:
    #   Internal Function called from na.omit.timSeries()

    # Arguments:

    # FUNCTION:

    # Extract matrix:
    x = object
    stopifnot (is.matrix(x))

    # Match Arguments:
    method <- match.arg(method)
    interp <- match.arg(interp)

    # Handle NAs:
    if (method == "r") {
        # Remove NAs:
        x <- na.omit(x)
    } else if (method == "z") {
        # Substitute NAs by Zero's:
        x[is.na(x)] <- 0
    } else if (substr(method, 1, 1) == "i") {
        # Interpolate:
        interp = match.arg(interp)
        f = 0
        if (interp == "before") {
            interp = "constant"
            f = 0
        }
        if (interp == "after") {
            interp = "constant"
            f = 1
        }
        n = nrow(x)
        for (i in 1:ncol(x)) {
            y = x[, i]
            idy = (1:n)[!is.na(y)]
            y = approx(idy, y[idy], 1:n, method = interp, f = f)$y
            x[, i] = y
        }
        modID = FALSE
        if (method == "ir") {
            # Remove Start and End NAs:
            x = na.omit(x)
        } else if (method == "iz") {
            # Set Start and End NAs to Zero:
            x[is.na(x)] = 0
        } else if (method == "ie") {
            n = nrow(x)
            for (k in 1:ncol(x)) {
                y  = x[, k]
                if (NA %in% y) {
                    start = sum(cumprod(is.na(y)))
                    if (start > 0) for (i in start:1) y[i] = y[i+1]
                    end = n+1 - sum(cumprod(rev(is.na(y))))
                    if (end <= n) for (i in end:n) y[i] = y[i-1]
                    x[, k] = y
                }
            }
        }
    }

    # Add Control:
    if (substr(method, 1, 1) == "i") {
        attr(x, "control") = c(method = method, interp = interp)
    } else {
        attr(x, "control") = c(method = method)
    }

    # Return Value:
    x
}


################################################################################


removeNA <- 
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Remove NA's from objects which can be transformed to a matrix

    # Arguments:
    #   x - an object which can be transformed to a matrix

    # FUNCTION:

    na.omit(x, ...)
}


# ------------------------------------------------------------------------------


substituteNA <- 
    function(x, type = c("zeros", "mean", "median"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Imputes missing data by zeros, the median or the
    #   mean values of all matrix elements

    # Arguments:
    #   x - an object which can be transformed to a matrix
    #   type - method specifies the substitution method to be
    #       used. Choices are "zeros", "mean", or "constant"

    # FUNCTION:
    if (!inherits(x, "matrix") && !inherits(x, "timeSeries"))
        x <- as(x, "matrix")

    # Type:
    type <- match.arg(type)
    ans <- switch(type,
        "zeros" = apply(x, 2, function(z) {
            z[is.na(z)] <- 0
            z}),
        "median" = apply(x, 2, function(z) {
            z[is.na(z)] = median(z, na.rm = TRUE)
            z}),
        "mean" = apply(x, 2, function(z) {
            z[is.na(z)] = mean(z, na.rm = TRUE)
            z}))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


interpNA <- 
    function(x, method = c("linear", "before", "after"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Interpolates missing values in a matrix object

    # Arguments:
    #   x - a numeric vector or time series object of class 'ts'.
    #   method - the method how to interpolate the vector, one of
    #       the applied vector strings: "linear", "before" or
    #       after.

    # Details:
    #   To interpolate the function 'approx' is used.

    # Value:
    #   Returns a vector or time series object where the missing
    #   values are interpolated.

    # FUNCTION:
    if (!inherits(x, "matrix") && !inherits(x, "timeSeries"))
        x <- as(x, "matrix")

    # Internal Function:
    interpVectorNA <- function(x, method, f) {
        n <- length(x)
        idx <- (1:n)[!is.na(x)]
        x <- approx(idx, x[idx], 1:n, method = method, f = f)$y
        x}

    # Select Method:
    method = method[1];
    f = 0
    if (method == "before") {
        method = "constant"
        f = 0
    }
    if (method == "after") {
        method = "constant"
        f = 1
    }

    # For each Column:
    for (i in 1:ncol(x)) {
        x[, i] = interpVectorNA(x[, i], method, f)
    }

    # Return Value:
    x
}


################################################################################

