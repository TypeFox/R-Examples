
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
# FUNCTION:                DESCRIPTION:
#   Zone                            Offset   Where
#   ----                            ------   ------
#    EET     Eastern European        +2      Finland, Eastern Europe
#    CET     Central European        +1      Western Europe, Sweden
#    GMT     Greenwich Mean          none    United Kingdom, Portugal
#    AST     Atlantic Standard       -4      Halifax
#    EST     Eastern Standard        -5      NY, DC, Toronto, Montreal
#    CST     Central Standard        -6      Chicago, Houston, Winnipeg
#    MST     Mountain Standard       -7      Denver, Calgary, Edmonton
#    PST     Pacific Standard        -8      LA, San Fransisco, Vancouver
#
# source http://www.stacken.kth.se/~kvickers/timezone.html
#
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
EET <- function()
    data.frame(EET = "1902-01-01 00:00:00",
               offSet = 2 * 3600,
               isdst = 0,
               TimeZone = "EET",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CET <- function()
    data.frame(CET = "1902-01-01 00:00:00",
               offSet = 1 * 3600,
               isdst = 0,
               TimeZone = "CET",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
AST <- function()
    data.frame(AST = "1902-01-01 00:00:00",
               offSet = -4 * 3600,
               isdst = 0,
               TimeZone = "AST",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
EST <- function()
    data.frame(EST = "1902-01-01 00:00:00",
               offSet = -5 * 3600,
               isdst = 0,
               TimeZone = "EST",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CST <- function()
    data.frame(CST = "1902-01-01 00:00:00",
               offSet = -6 * 3600,
               isdst = 0,
               TimeZone = "CST",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
MST <- function()
    data.frame(MST = "1902-01-01 00:00:00",
               offSet = -7 * 3600,
               isdst = 0,
               TimeZone = "MST",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
PST <- function()
    data.frame(PST = "1902-01-01 00:00:00",
               offSet = -8 * 3600,
               isdst = 0,
               TimeZone = "PST",
               numeric = as.numeric(as.POSIXct("1902-01-01 00:00:00")),
               stringsAsFactors = FALSE)
