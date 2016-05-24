#######################################################################
# rEMM - Extensible Markov Model (EMM) for Data Stream Clustering in R
# Copyrigth (C) 2011 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## helper
.installed <- function(pkg) !is(try(installed.packages()[pkg,],
        silent=TRUE), "try-error")


.get_parameters <- function(p, parameter) {
    if(!is.null(parameter) && length(parameter) != 0) {
        o <- pmatch(names(parameter), names(p))

        if(any(is.na(o)))
        stop(sprintf(ngettext(length(is.na(o)),
                    "Unknown option: %s",
                    "Unknown options: %s"),
                paste(names(parameter)[is.na(o)],
                    collapse = " ")))

        p[o] <- parameter
    }

    p
}


