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
# FUNCTION:                  DESCRIPTION:
#  na.contiguous,timeSeries  Finds the longest consecutive of non-missing values
################################################################################


setMethod("na.contiguous", "timeSeries", 
    function(object, ...)
{
    # A function imlemented by Diethelm Wuertz and Yohan Chalabi
      
    # Description:
    #   Finds the longest consecutive of non-missing values
    
    # Details:
    #     adapted stats:::na.contingous.default to timeSeries objects
    #     Yohan Chalabi
    
    # FUNCTION:
    
    good <- apply(!is.na(object), 1L, all)
    if (!sum(good))
        stop("all times contain an NA")
    tt <- cumsum(!good)
    ln <- sapply(0:max(tt), function(i) sum(tt == i))
    seg <- (seq_along(ln)[ln == max(ln)])[1L] - 1
    keep <- (tt == seg)
    st <- min(which(keep))
    if (!good[st])
        st <- st + 1
    en <- max(which(keep))
    omit <- integer(0L)
    n <- NROW(object)
    if (st > 1)
        omit <- c(omit, 1L:(st - 1))
    if (en < n)
        omit <- c(omit, (en + 1):n)
    if (length(omit)) {
        object <- object[st:en, ]
        attr(omit, "class") <- "omit"
        attr(object, "na.action") <- omit
    }
    
    # Return Value:
    object
})


################################################################################

