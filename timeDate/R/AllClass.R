
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

# Copyrights (C)
# for this R-port:
#   1999 - Diethelm Wuertz, GPL
#   2007 - Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                GENERATION OF TIMEDATE OBJECTS:
#  setClass                 'timeDate' S4 Class representation
#  setMethod inti           'initialize', 'timeDate
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setClass("timeDate",
         # A class implemented by Diethelm Wuertz and Yohan Chalabi

         # Description:
         #   Class representatation for 'timeDate' Objects.

         # CLASS:
         representation(Data = "POSIXct",
                        format = "character",
                        FinCenter = "character"),
         validity = function(object) {
             if(!identical(attr(object@Data, "tzone"), "GMT"))
                 return("@Data must be in \"GMT\" timezone.")
             if(!is.numeric(unclass(object@Data)))
                 return("unclass(@Data) should be of class \"numeric\".")
             ## else TRUE
             TRUE
         })

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("initialize", "timeDate", function(.Object, ...) {

    .Object <- callNextMethod()

    # if not arguments are passed in ..., do not try to define format
    # of @Data
    if (length(list(...))) {

        # ISO Date/Time Format:
        isoDate   <- "%Y-%m-%d"
        isoFormat <- "%Y-%m-%d %H:%M:%S"

        # extract numerical value
        num <- c(unclass(.Object@Data))

        if (all(is.na(num))) {
            # no need to look for a format if @Data has only NA's
            .Object@format <- character(1)
        } else {
            # convert - DST
            num <- .formatFinCenterNum(num, .Object@FinCenter, "gmt2any")

            # check if num is a multiple of days
            test <- !(abs(num %% 86400) > 0)
            .Object@format <- ifelse(all(na.omit(test)), isoDate, isoFormat)
        }
    }
    .Object
})

################################################################################
