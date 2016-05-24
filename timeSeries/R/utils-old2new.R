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


.old2newTimeSeries <- 
function(x)
{

    # Version 1:
    if ("Data" %in% slotNames(x)) {
        data <- x@Data
        charvec <- timeDate(x@positions, zone = x@FinCenter,
                            FinCenter = x@FinCenter)
        units <- x@units
        recordIDs <- x@recordIDs
        title <- x@title
        documentation <- x@documentation
        x <- timeSeries(data = data,
                        charvec = charvec,
                        units = units,
                        recordIDs = recordIDs,
                        title = title,
                        documentation = documentation)
    }

    # Version 2:
    if ((".Data" %in% slotNames(x)) && is.character(x@positions)) {
        data <- x@.Data
        charvec <- timeDate(x@positions, zone = x@FinCenter,
                            FinCenter = x@FinCenter)
        units <- x@units
        recordIDs <- x@recordIDs
        title <- x@title
        documentation <- x@documentation
        x <- timeSeries(data = data,
                        charvec = charvec,
                        units = units,
                        recordIDs = recordIDs,
                        title = title,
                        documentation = documentation)
    }

    x
}


# ------------------------------------------------------------------------------


## # Example
## library(timeSeries)
## setwd("~/r/fPortfolio/data")
## rda <- dir()
## sapply(rda, .old2newRda, suffix = "")


.old2newRda <- function(file, suffix = "_new")
{
    stopifnot(length(file) == 1)
    local({
        load(file)
        nm <- ls()
        lold <- mget(nm, envir = environment(NULL))
        test <- sapply(lold, is.timeSeries)
        lold <- lold[test]
        lnew <- lapply(lold, .old2newTimeSeries)
        objects <- names(lold)
        for (nm in objects)
            assign(nm, lnew[[nm]])
        newFile <- paste(file, suffix, sep = "")
        save(list = objects, file = newFile)
    })

    invisible(TRUE)
}


################################################################################

