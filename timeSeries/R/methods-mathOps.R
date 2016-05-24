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
#  Ops,timeSeries            Returns group 'Ops' for a 'timeSeries' object
#  cummax                    Returns cumulated maxima
#  cummin                    Returns cumulated minima
#  cumprod                   Returns cumulated products
#  cumsum                    Returns cumulated sums
#  ##diff,timeSeries           Differences a timeSeries object
#  ##scale,timeSeries          Scales a timeSeries object
#  quantile,timeSeries       Samples qunatiles of a timeSeries object
################################################################################


setMethod("Ops", c("vector", "timeSeries"),
    function(e1, e2)
    {
          lattrs <- attributes(e2)
          e2 <- getDataPart(e2)
          value <- callGeneric(e1, e2)
          if (identical(dim(value), dim(e2))) {
              attributes(value) <- lattrs
              value <- asS4(value, TRUE)
          }
          value
    }
)


# ------------------------------------------------------------------------------


setMethod("Ops", c("array", "timeSeries"),
    function(e1, e2)
    {
          e1 <- as.vector(e1)
          lattrs <- attributes(e2)
          e2 <- getDataPart(e2)
          value <- callGeneric(e1, e2)
          if (identical(dim(value), dim(e2))) {
              attributes(value) <- lattrs
              value <- asS4(value, TRUE)
          }
          value
    }
)


# ------------------------------------------------------------------------------


setMethod("Ops", c("ts", "timeSeries"),
    function(e1, e2)
    {
          e1 <- as(e1, "matrix")
          lattrs <- attributes(e2)
          e2 <- getDataPart(e2)
          value <- callGeneric(e1, e2)
          if (identical(dim(value), dim(e2))) {
              attributes(value) <- lattrs
              value <- asS4(value, TRUE)
          }
          value
    }
)


# ------------------------------------------------------------------------------


setMethod("Ops", c("timeSeries", "vector"),
    function(e1, e2)
    {
          lattrs <- attributes(e1)
          e1 <- getDataPart(e1)
          value <- callGeneric(e1, e2)
          if (identical(dim(value), dim(e1))) {
              attributes(value) <- lattrs
              value <- asS4(value, TRUE)
          }
          value
    }
)


# ------------------------------------------------------------------------------


setMethod("Ops", c("timeSeries", "array"),
    function(e1, e2)
    {
          lattrs <- attributes(e1)
          e1 <- getDataPart(e1)
          value <- callGeneric(e1, e2)
          if (identical(dim(value), dim(e1))) {
              attributes(value) <- lattrs
              value <- asS4(value, TRUE)
          }
          value
    }
)


# ------------------------------------------------------------------------------


setMethod("Ops", c("timeSeries", "ts"),
    function(e1, e2)
    {
          lattrs <- attributes(e1)
          e1 <- getDataPart(e1)
          e2 <- as(e2, "matrix")
          value <- callGeneric(e1, e2)
          if (identical(dim(value), dim(e1))) {
              attributes(value) <- lattrs
              value <- asS4(value, TRUE)
          }
          value
    }
)


# ------------------------------------------------------------------------------


setMethod("Ops", c("timeSeries", "timeSeries"),
    function(e1, e2)
    {
        # Note keep recordIDs of e1 only
        # check if conformable arrays
        if (!identical(dim(e1), dim(e2)))
            stop("non-conformable arrays")
        # check if positions are identical
        if (!identical(e1@positions, e2@positions))
            stop("positions slot do not match")
        lattrs <- attributes(e1)
        e1 <- getDataPart(e1)
        e2 <- getDataPart(e2)
        value <- callGeneric(e1, e2)
        if (identical(dim(value), dim(e1))) {
            attributes(value) <- lattrs
            value <- asS4(value, TRUE)
        }
        value
    }
)


# ------------------------------------------------------------------------------


setMethod("cummax", "timeSeries",
    function(x) callGeneric(getDataPart(x)))


# ------------------------------------------------------------------------------


setMethod("cummin", "timeSeries",
    function(x) callGeneric(getDataPart(x)))


# ------------------------------------------------------------------------------


setMethod("cumprod", "timeSeries",
    function(x) callGeneric(getDataPart(x)))


# ------------------------------------------------------------------------------


setMethod("cumsum", "timeSeries",
    function(x) callGeneric(getDataPart(x)))


# ------------------------------------------------------------------------------


## setMethod("diff", "timeSeries", function(x, ...) {
##    x <- getDataPart(x)
##    callGeneric()
## })


# ------------------------------------------------------------------------------


## setMethod("scale", "timeSeries", function(x, ...) {
##    x <- getDataPart(x)
##    callGeneric()
## })


# ------------------------------------------------------------------------------


setMethod("quantile", "timeSeries", function(x, ...) {
    x <- getDataPart(x)
    callGeneric()
})


################################################################################
