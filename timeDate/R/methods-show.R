
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  show.timeDate             Prints 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("show", "timeDate", function (object)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # when creating empty new("timeDate")
    if (!length(slot(object, "Data")))
        return(str(object))

    # Check records to get printed:
    maxRmetrics <- as.numeric(getRmetricsOptions("max.print"))
    maxR <- as.numeric(getOption("max.print"))
    max <- min(na.omit(c(maxRmetrics, maxR, Inf)))
    #-> Inf to cast case when maxRmetrics and maxR are NULL

    if (ptest <- ((omitted <- length(object) - max) > 0))
        object <- object[seq.int(max)]

    output <- format(object)
    layout <- paste0("[", output, "]")
    names(layout) <- names(output)

    # Print Results:
    cat(finCenter(object), "\n", sep = "")
    print(layout, quote = FALSE)

    # print message
    if (ptest)
        cat(gettextf("...\n [ reached getRmetricsOption('max.print') | getOption('max.print') -- omitted %i rows ]]\n", omitted))

    # Return Value:
    invisible(NULL) # 'show' returns an invisible 'NULL'. (cf. ?show)
})


################################################################################
