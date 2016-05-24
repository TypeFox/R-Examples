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
#  model.frame.default       Allows to use model.frame for "timeSeries"
################################################################################


# YC : remove model.frame because more problems than benefits. Rely on
# default model.frame as long as as.data.frame.timeSeries works in
# 'base' function model.frame.default


## setMethod("model.frame.default", signature(data = "timeSeries"),
##           function(formula, data = NULL,
##                    subset = NULL, na.action = na.fail,
##                    drop.unused.levels = FALSE, xlev = NULL, ...)
##       {
##           # A function implemented by Diethelm Wuertz

##           # Description:
##           #   Extracting the Environment of a Model Formula

##           # Arguments:
##           #   formula - a model formula
##           #   data - a 'timeSeries' object

##           # Details:
##           #   Allows to use model.frame() for "timeSeries" objects.

##           # Examples:
##           #   x = as.timeSeries(data(msft.dat))[1:12, ]
##           #   model.frame( ~ High + Low, data = x)
##           #   model.frame(Open ~ High + log(Low), data = x)

##           # FUNCTION:
##           data <- as(data, "data.frame")

## ###           model.frame.default(formula, data,
## ###                               subset, na.action,
## ###                               drop.unused.levels,
## ###                               xlev, ...)

##           model.frame(formula, data, ...)
##       })

## ## ###           # Create Model Frame:
## ## ###           format <- data@format
## ## ###           FinCenter <- finCenter(data)
## ## ###           recordIDs <- data@recordIDs
## ## ###           title <- data@title

## ##           data <- as(data, "data.frame")
## ##           Model <- model.frame(formula, data, ...)
## ##           #-> should be in parent.frame?

## ## ###           recordIDs <-
## ## ###               if (NROW(Model) == NROW(recordIDs))
## ## ###                   recordIDs
## ## ###               else
## ## ###                   data.frame()

## ## ###           # Convert to timeSeries:
## ## ###           ans <- timeSeries(data = as.matrix(Model),
## ## ###                             charvec = rownames(Model),
## ## ###                             units = colnames(Model),
## ## ###                             format = format,
## ## ###                             FinCenter = FinCenter,
## ## ###                             recordIDs = recordIDs,
## ## ###                             title = title,
## ## ###                             documentation = description()
## ## ###                             )

## ## ###           # Return value:
## ## ###           ans
## ##           Model
## ##       })


################################################################################

