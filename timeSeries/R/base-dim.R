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
# S4 METHOD:                DIM OPERATIONS ON DATA:
#  dim,timeSeries            Returns dimension of a 'timeSeries' object
#  dim<-,timeSeries          Assigns dimension of a 'timeSeries' object
#  dimnames,timeDSeries      Returns dimension names of a 'timeSeries' object
#  dimnames<-,timeSeries     Assign dimension names of a 'timeSeries' object
#  colnames,timeSeries       Return column names to a 'timeSeries' object
#  rownames,timeSeries       Return row names to a 'timeSeries' object
#  colnames<-,timeSeries     Assigns column names to a 'timeSeries' object
#  rownames<-,timeSeries     Assigns row names to a 'timeSeries' object
#  names,timeSeries          Return column names of a 'timeSeries' object
#  names<.,timeSeries        Assigns column names of a 'timeSeries' object
################################################################################


# Base Functions:

    # Generate from Matrix:
    # edhec.tS = timeSeries(edhec.mat, rownames(edhec.mat))
    # edhec.ts = ts(edhec.mat, start = c(1997, 1), frequency = 12)

    # Univariate time Series:
    # edhec1.tS = edhec.tS[, 1]

    #   dim
    #                       dim(edhec.tS)                       # 20 4
    #                       dim(edhec1.tS)                      # 20 1


    #   DIM
    #                       DIM = function(x) {c(NROW(x), NCOL(x))}
    #                       DIM(edhec.tS)                       # 20 4
    #                       DIM(edhec1.tS)                      # 20 1


    #   length
    #                       length(edhec.tS)                    # 1

    #
    #   LENGTH
    #                       LENGTH = function(x) NROW(x)
    #                       LENGTH(edhec.tS)                    # 20
    #                       LENGTH(edhec1.tS)                   # 20


    #
    #   ncol / nrow
    #                       ncol(edhec.tS)                      # 4

    #
    #                       ncol(edhec1.tS)                     # 1

    #
    #  NCOL / NRWO
    #                       NCOL(edhec.tS)                      # 4

    #
    #                       NCOL(edhec1.tS)                     # 1

    #
    #  isUnivariate
    #                       isUnivariate = function(x) NCOL(x) == 1
    #                       isUnivariate(edhec.tS)
    #                       isUnivariate(edhec1.tS)


    #
    # isMultivariate        # Just Negation of isUnivariate
    #
    #
    #

# ------------------------------------------------------------------------------


# length
# dim
# ncol
# nrow


# LENGTH
# DIM
# NCOL
# NROW

# ------------------------------------------------------------------------------


# Note it is faster to access attribute rather than accessing @.Data


setMethod("dim", "timeSeries", function(x) attr(x, "dim"))

# This should make functions like
# model.response(model.frame(dummySeries() ~1)) work


setReplaceMethod("dim", "timeSeries",
    function(x, value)
    {
        # dim(x) <- NULL returns a vector
        if (is.null(value))
            return(as.vector(x))
        else
            x #<< returns same object : # setting new dim
            # is forbidden for a timeSeries object
    }
)


# ------------------------------------------------------------------------------


# colnames - faster to have dedicated method than relying on dimnames[[2]]


setMethod("colnames", "timeSeries", # "signalSeries",
    function(x, do.NULL = TRUE, prefix = "col") x@units
)


# ------------------------------------------------------------------------------


# rownames 


## setMethod("rownames", "signalSeries",
##           function (x, do.NULL = TRUE, prefix = "row") NULL)

## setMethod("rownames", "timeSeries",
##           function (x, do.NULL = TRUE, prefix = "row") as.character(time(x)))


setMethod("rownames", "timeSeries",
    function (x, do.NULL = TRUE, prefix = "row")
    {
        if (length(x@positions) > 0)
            as.character(time(x))
        else
            NULL
    }
)


# ------------------------------------------------------------------------------


setMethod("dimnames", "timeSeries", # "signalSeries",
    function(x)
    {
        list(rownames(x),colnames(x))
    }
)


# ------------------------------------------------------------------------------


setMethod("colnames<-", "timeSeries",
    function(x, value)
    {
        units <- as.character(value)

        if(!length(units))
            if (x@format == "counts")
                  units <- paste("SS", seq(NCOL(x)), sep = ".")
              else
                  units <- paste("TS", seq(NCOL(x)), sep = ".")

        if (length(units) != NCOL(x))
            stop("length of 'colnames' not equal to array extent",call.=FALSE)

        x@units <- units
        colnames(x@.Data) <- units

        x
    }
)


# ------------------------------------------------------------------------------


setMethod("rownames<-", c("timeSeries", "timeDate"), #c("signalSeries", "timeDate"),
    function (x, value)
    {
        .timeSeries(
            data = getDataPart(x),
            charvec = as.numeric(value, "sec"),
            units = colnames(x),
            format = value@format,
            FinCenter = value@FinCenter,
            recordIDs = x@recordIDs,
            title = x@title,
            documentation = x@documentation)
    }
)


# ------------------------------------------------------------------------------


setMethod("rownames<-", "timeSeries", # "signalSeries",
    function (x, value)
    {
        # if charvec NULL returns a signal series
        if (is.null(value))
            return(.signalSeries(data = getDataPart(x),
                                 units = colnames(x),
                                 recordIDs = x@recordIDs,
                                 title = x@title,
                                 documentation = x@documentation))

        # coerce charvec to timeDate
        charvec <- timeDate(charvec = value)

        if (any(is.na(charvec)))
            # Note : there is already a warning in timeDate if there are NA's
            .signalSeries(data = getDataPart(x),
                          units = colnames(x),
                          recordIDs = x@recordIDs,
                          title = x@title,
                          documentation = x@documentation)
        else
            .timeSeries(data = getDataPart(x),
                        charvec = as.numeric(charvec, "sec"),
                        units = colnames(x),
                        format = charvec@format,
                        FinCenter = charvec@FinCenter,
                        recordIDs = x@recordIDs,
                        title = x@title,
                        documentation = x@documentation)
    }
)


# ------------------------------------------------------------------------------


setMethod("dimnames<-", c("timeSeries", "list"), # c("signalSeries", "list"),
    function(x, value)
    {
        rownames(x) <- value[[1]]
        colnames(x) <- value[[2]]

        x
    }
)


# ------------------------------------------------------------------------------
# important for completion with $


setMethod("names", "timeSeries", # "signalSeries",
    function(x) c(colnames(x), names(x@recordIDs)))


setReplaceMethod("names", "timeSeries", # "signalSeries",
function(x, value) {
   nc <- ncol(x)
   nv <- length(value)
   nr <- length(x@recordIDs)

   # Note that using [][] ensure that length of the
   # names are equal to array extent
   colnames(x) <- value[seq.int(nv)][seq.int(nc)]
   if (nv > nc)
       names(x@recordIDs) <- value[-seq.int(nc)][seq.int(nr)]
   x
})


################################################################################
